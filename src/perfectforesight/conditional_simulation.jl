function set_future_information!(y::Vector{Float64}, x::Vector{Float64}, context, periods, infoperiod)
    m = context.models[1]
    endogenous_nbr = m.endogenous_nbr
    exogenous_nbr = m.exogenous_nbr
    flips = context.work.scenario[infoperiod]
    symboltable = context.symboltable
    
    for (period, f1) in flips
        p = length(infoperiod:period)
        for (s, f2) in f1
            ss = string(s)
            if is_exogenous(ss, symboltable)
                idx = (p - 1)*exogenous_nbr + symboltable[ss].orderintype
                x[idx] = f2[1]
            elseif is_endogenous(ss, symboltable)
                idy = (p - 1)*endogenous_nbr + symboltable[ss].orderintype
                y[idy] = f2[1]
            end
        end
    end
end
        
function perfectforesight_core_conditional!(
    perfect_foresight_ws::PerfectForesightWs,
    context::Context,
    periods::Int64,
    y0::AbstractVector{Float64},
    initialvalues::Vector{<:Real},
    terminalvalues::Vector{<:Real},
    linear_solve_algo::LinearSolveAlgo,
    dynamic_ws::DynamicWs,
    flipinfo::FlipInformation,
    infoperiod
)
    m = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    residuals = zeros(periods * m.endogenous_nbr)
    dynamic_variables = dynamic_ws.dynamic_variables
    temp_vec = dynamic_ws.temporary_values
    steadystate = results.trends.endogenous_steady_state
    params = work.params
    JJ = perfect_foresight_ws.J

    exogenous = perfect_foresight_ws.x
    ws_threaded = [
        Dynare.DynamicWs(
            m.endogenous_nbr,
            m.exogenous_nbr,
            length(temp_vec),
            m.dynamic_g1_sparse_colptr,
            m.dynamic_g1_sparse_rowval
        ) for i = 1:Threads.nthreads()
    ]
    nzval = similar(m.dynamic_g1_sparse_rowval, Float64)
    nzval1 = similar(nzval)
    
    f! = make_pf_residuals(
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        m,
        periods,
        temp_vec,
        perfect_foresight_ws.permutationsR,
        flipinfo.ix_stack
    )

    J! = make_pf_jacobian(
        DFunctions.dynamic_derivatives!,
        initialvalues,
        terminalvalues,
        dynamic_variables,
        exogenous,
        periods,
        temp_vec,
        params,
        steadystate,
        m.dynamic_g1_sparse_colptr,
        nzval,
        m.endogenous_nbr,
        m.exogenous_nbr,
        perfect_foresight_ws.permutationsJ,
        flipinfo,
        nzval1
    )

    df = OnceDifferentiable(f!, J!, y0, residuals, JJ)
    @debug "$(now()): start nlsolve"

    set_future_information!(y0, exogenous, context, periods, infoperiod)
    flip!(y0, exogenous, flipinfo.ix_stack)
    f!(residuals, y0)
    if linear_solve_algo == pardiso
        @show "Pardiso"
        if isnothing(Base.get_extension(Dynare, :PardisoSolver))
            error("You must load Pardiso with 'using MKL, Pardiso'")
        end
        ls1!(x, A, b) = linear_solver!(PardisoLS(), x, A, b)
        res = nlsolve(df, y0, method = :robust_trust_region, show_trace = true, ftol=cbrt(eps()), linsolve = ls1!)    
    else
        ls2!(x, A, b) = linear_solver!(IluLS(), x, A, b)
        res = nlsolve(df, y0, method = :robust_trust_region, show_trace = false, ftol=cbrt(eps()), linsolve = ls2!)
    end
    print_nlsolver_results(res)
    @debug "$(now()): end nlsolve"
    flip!(res.zero, exogenous, flipinfo.ix_stack)
    make_simulation_results!(context::Context, res.zero, exogenous, terminalvalues, periods)
end

function make_pf_residuals(
    initialvalues::AbstractVector{T},
    terminalvalues::AbstractVector{T},
    exogenous::AbstractVector{T},
    dynamic_variables::AbstractVector{T},
    steadystate::AbstractVector{T},
    params::AbstractVector{T},
    m::Model,
    periods::Int,
    temp_vec::AbstractVector{T},
    permutations::Vector{Tuple{Int64,Int64}},
    ix_stack::SparseVector{Int, Int}
        ) where T <: Real
    function f!(residuals::AbstractVector{T}, y::AbstractVector{T})
        flip!(y, exogenous, ix_stack)
        get_residuals!(
            residuals,
            vec(y),
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            permutations = permutations
        )
        flip!(y, exogenous, ix_stack)
        return residuals
    end
    return f!
end

function flip!(y, x, ix_stack::SparseVector{Int, Int})
    I, J = findnz(ix_stack)
    for (iy, ix) in zip(I, J)
        y[iy], x[ix] = x[ix], y[iy]
    end
end

function make_pf_jacobian(
    dynamic_derivatives!::Function,
    initialvalues::AbstractVector{T},
    terminalvalues::AbstractVector{T},
    dynamic_variables::AbstractVector{T},
    exogenous::AbstractVector{T},
    periods::N,
    temp_vec::AbstractVector{T},
    params::AbstractVector{T},
    steadystate::AbstractVector{T},
    dynamic_g1_sparse_colptr::AbstractVector{N},
    nzval::AbstractVector{T},
    endogenous_nbr::N,
    exogenous_nbr::N,
    permutations,
    flipinfo::FlipInformation,
    nzval1::AbstractVector{T},
) where {T <: Real, N <: Integer}
    function J!(
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        y::AbstractVector{Float64},
    )
        updateJacobian!(
            A,
            dynamic_derivatives!,
            y,
            initialvalues,
            terminalvalues,
            dynamic_variables,
            exogenous,
            periods,
            temp_vec,
            params,
            steadystate,
            dynamic_g1_sparse_colptr,
            nzval,
            endogenous_nbr,
            exogenous_nbr,
            permutations,
            flipinfo,
            nzval1
        )
        return A
    end
end

function adjust_colptr_size(iy_period::SparseVector{Int, Int}, colptr::Vector{Int}, endogenous_nbr::Int, periods::Int)
    size_adjustment = 0
    I, J = findnz(iy_period)
    for (iy2, iy1) in zip(I, J)
        if iy2 <= endogenous_nbr
            size_adjustment += (1 - (colptr[iy1 + 1] - colptr[iy1])
                                - (colptr[iy1 + endogenous_nbr + 1] - colptr[iy1 + endogenous_nbr]))

        elseif iy2 > (periods - 1)*endogenous_nbr
            size_adjustment += (1 - (colptr[iy1 + endogenous_nbr + 1] - colptr[iy1 + endogenous_nbr])
                                - (colptr[iy1 + 2*endogenous_nbr + 1] - colptr[iy1 + 2*endogenous_nbr]))
        else
            size_adjustment += (1 - (colptr[iy1 + 1] - colptr[iy1])
                                - (colptr[iy1 + endogenous_nbr + 1] - colptr[iy1 + endogenous_nbr])
                                - (colptr[iy1 + 2*endogenous_nbr + 1] - colptr[iy1 + 2*endogenous_nbr]))
        end
    end
    return size_adjustment
end

function makeJacobian(colptr, rowval, endogenous_nbr, periods, permutations,
    flipinfo::FlipInformation
    )

    flipped_variables = flipinfo.flipped_variables
    ix_period = flipinfo.ix_period
    iy_period = flipinfo.iy_period
    
    startv = colptr[endogenous_nbr + 1]
    endv = colptr[2*endogenous_nbr + 1]

    rowval_length = colptr[3*endogenous_nbr + 1] - 1
    nnz = (periods - 1)*rowval_length - startv + endv
    nnz += adjust_colptr_size(iy_period, colptr, endogenous_nbr, periods)
    # update nnz
    bigcolptr = Vector{Int64}(undef, periods*endogenous_nbr + 1)
    bigrowval = Vector{Int64}(undef, nnz)
    bignzval = similar(bigrowval, Float64)
    permutations1 = []

    # compute permutations for one period Jacobian
    if !isempty(permutations)
        (permutations1, rowval) = permutation(permutations, colptr, rowval)
    end

    r = 1
    c = 1
    # first periods
    bigcolptr[1] = 1
    for i in 1:endogenous_nbr
        c += 1
        offset = 0
        if flipped_variables[c - 1]
            ix = ix_period[c - 1]
            bigcolptr[c] = bigcolptr[c - 1] + colptr[ix + 1] - colptr[ix]
            r = bigindex!(bigrowval, 
                          r,
                          rowval,
                          colptr[ix],
                          colptr[ix + 1] - 1,
                          offset)
        else
            bigcolptr[c] = (bigcolptr[c - 1]
                            + colptr[endogenous_nbr + i + 1]
                            - colptr[endogenous_nbr + i]
                            + colptr[i + 1]
                            - colptr[i])
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[endogenous_nbr + i],
                         colptr[endogenous_nbr + i + 1] - 1,
                         offset)
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[i],
                         colptr[i + 1] - 1,
                         offset + endogenous_nbr)
        end
    end

    # intermediary periods
    for p in 2:periods - 1
        for i in 1:endogenous_nbr
            c += 1
            offset = (p - 2)*endogenous_nbr
            # duplicate for flipped variable
            if flipped_variables[c - 1]
                ix = ix_period[c - 1]
                bigcolptr[c] = bigcolptr[c - 1] + colptr[ix + 1] -colptr[ix]
                r = bigindex!(bigrowval, 
                              r,
                              rowval,
                              colptr[ix],
                              colptr[ix + 1] - 1,
                              offset + endogenous_nbr)
            else
                bigcolptr[c] = (bigcolptr[c - 1]
                                + colptr[2*endogenous_nbr + i + 1]
                                - colptr[2*endogenous_nbr + i]
                                + colptr[endogenous_nbr + i + 1]
                                - colptr[endogenous_nbr + i]
                                + colptr[i + 1]
                                - colptr[i])
                r= bigindex!(bigrowval,
                             r,
                             rowval,
                             colptr[2*endogenous_nbr + i],
                             colptr[2*endogenous_nbr + i + 1] - 1,
                             offset)
                r= bigindex!(bigrowval,
                             r,
                             rowval,
                             colptr[endogenous_nbr + i],
                             colptr[endogenous_nbr + i + 1] - 1,
                             offset + endogenous_nbr)
                r= bigindex!(bigrowval,
                             r,
                             rowval,
                             colptr[i],
                             colptr[i + 1] - 1,
                             offset + 2*endogenous_nbr)
            end
        end
        offset += endogenous_nbr
    end
    #terminal period
    for i in 1:endogenous_nbr
        c += 1
        offset = (periods - 2)*endogenous_nbr
        if flipped_variables[c - 1]
            ix = ix_period[c - 1]
            bigcolptr[c] = bigcolptr[c - 1] + colptr[ix + 1] -colptr[ix]
            r = bigindex!(bigrowval, 
                          r,
                          rowval,
                          colptr[ix],
                          colptr[ix + 1] - 1,
                          offset + endogenous_nbr)
        else
            bigcolptr[c] = (bigcolptr[c - 1]
                            + colptr[2*endogenous_nbr + i + 1]
                            - colptr[2*endogenous_nbr + i]
                            + colptr[endogenous_nbr + i + 1]
                            - colptr[endogenous_nbr + i])
            r1 = colptr[endogenous_nbr + i]
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[2*endogenous_nbr + i],
                         colptr[2*endogenous_nbr + i + 1] - 1,
                         offset)
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[endogenous_nbr + i],
                         colptr[endogenous_nbr + i + 1] - 1,
                         offset + endogenous_nbr)
        end
    end

    b1 = colptr[endogenous_nbr + 1]
    b2 = (periods - 1)*rowval_length + colptr[2*endogenous_nbr + 1] - 1
    nperiods = periods*endogenous_nbr
    return (SparseMatrixCSC(nperiods, nperiods, bigcolptr, bigrowval, bignzval), permutations1)
end

function updateJacobian!(J::SparseMatrixCSC,
                         G1!,
                         endogenous::AbstractVector{<: Real},
                         initialvalues::AbstractVector{<: Real},
                         terminalvalues::AbstractVector{<: Real},
                         dynamic_variables::AbstractVector{<: Real},
                         exogenous::AbstractVector{<: Real},
                         periods,
                         temporary_var::AbstractVector{<: Real},
                         params::AbstractVector{<: Real},
                         steady_state::AbstractVector{<: Real},
                         colptr::AbstractVector{Int},
                         nzval::AbstractVector{<: Real},
                         endogenous_nbr,
                         exogenous_nbr,
                         permutations,
                         flipinfo::FlipInformation,
                         ws::AbstractVector{<: Real})
    flipped_variables = flipinfo.flipped_variables
    ix_period = flipinfo.ix_period
    bigcolptr = J.colptr
    offset = 1
    rx = 1:exogenous_nbr
    @views begin
        copyto!(dynamic_variables, initialvalues)
        copyto!(dynamic_variables, endogenous_nbr + 1, endogenous, 1, 2*endogenous_nbr)
        G1!(temporary_var, nzval, dynamic_variables, exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        oy = endogenous_nbr
        ox = exogenous_nbr
        for c in 1:2*endogenous_nbr
            if flipped_variables[c] && c <= endogenous_nbr
                n = colptr[ix_period[c] + 1] - colptr[ix_period[c]]
                copyto!(J.nzval, bigcolptr[c], nzval, colptr[ix_period[c]], n)
            else
                k = bigcolptr[c]
                n = colptr[endogenous_nbr + c+1] - colptr[endogenous_nbr + c]
                copyto!(J.nzval, k, nzval, colptr[endogenous_nbr + c], n)
            end
        end
        ry = 1:3*endogenous_nbr
        rx = rx .+ ox
        G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        for c in 1:2*endogenous_nbr
            if flipped_variables[c] && c > endogenous_nbr
                n = colptr[ix_period[c] + 1] - colptr[ix_period[c]]
                copyto!(J.nzval, bigcolptr[c], nzval, colptr[ix_period[c]], n)
            else
                k = bigcolptr[c] + colptr[endogenous_nbr + c + 1] - colptr[endogenous_nbr + c]
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
        end
        for c in 2*endogenous_nbr + 1:3*endogenous_nbr
            k = bigcolptr[c]
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end

        for p in 1:periods - 3
            ry = ry .+ oy
            rx = rx .+ ox
            G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state)
            !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
            for c in 1:endogenous_nbr
                k = (bigcolptr[c + p*endogenous_nbr]
                     + colptr[c + 2*endogenous_nbr + 1] - colptr[c + 2*endogenous_nbr]
                     + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
            for c in endogenous_nbr + 1:2*endogenous_nbr
                c1 = c + p*endogenous_nbr
                if flipped_variables[c1]
                    n = colptr[ix_period[c1] + 1] - colptr[ix_period[c1]]
                    copyto!(J.nzval, bigcolptr[c1], nzval, colptr[ix_period[c1]], n)
                else
                    k = (bigcolptr[c1]
                         + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
                    n = colptr[c + 1] - colptr[c]
                    copyto!(J.nzval, k, nzval, colptr[c], n)
                end
            end
            for c in 2*endogenous_nbr + 1: 3*endogenous_nbr
                k = bigcolptr[c + p*endogenous_nbr]
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
        end

        rx = rx .+ ox
        copyto!(dynamic_variables, 1, endogenous, (periods - 2)*endogenous_nbr + 1, 2*endogenous_nbr)
        copyto!(dynamic_variables, 2*endogenous_nbr + 1, terminalvalues)
        G1!(temporary_var, nzval, dynamic_variables, exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        for c in 1:endogenous_nbr
            k = (bigcolptr[c + (periods - 2)*endogenous_nbr]
                 + colptr[c + 2*endogenous_nbr + 1] - colptr[c + 2*endogenous_nbr]
                 + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end
        for c in endogenous_nbr + 1:2*endogenous_nbr
            c1 = c + (periods - 2)*endogenous_nbr
            if flipped_variables[c1]
                n = colptr[ix_period[c1] + 1] - colptr[ix_period[c1]]
                copyto!(J.nzval, bigcolptr[c1], nzval, colptr[ix_period[c1]], n)
            else
                k = (bigcolptr[c1]
                     + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
        end
    end
    return J
end

