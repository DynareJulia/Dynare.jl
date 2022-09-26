using LinearAlgebra
using Suppressor

export get_de, get_abc, inverse_order_of_dynare_decision_rule

struct DynamicWs
    dynamic_variables::Vector{Float64}
    dynamic_variables2::Vector{Float64}
    exogenous_variables::Vector{Float64}
    residuals::Vector{Float64}
    derivatives::Vector{Vector{Float64}}
    temporary_values::Vector{Float64}
    nzval::Vector{Float64}
    function DynamicWs(
        endogenous_nbr::Int64,
        exogenous_nbr::Int64,
        dynamic_nbr::Int64,
        tmp_nbr::Int64,
        nzval_nbr::Int64
    )
        dynamic_variables = Vector{Float64}(undef, dynamic_nbr)
        dynamic_variables2 = Vector{Float64}(undef, 3*endogenous_nbr)
        exogenous_variables = Vector{Float64}(undef, exogenous_nbr)
        residuals = Vector{Float64}(undef, endogenous_nbr)
        derivatives = Vector{Vector{Float64}}(undef, 1)
        derivatives[1] = Vector{Float64}(undef, 0)
        temporary_values = Vector{Float64}(undef, tmp_nbr)
        nzval = Vector{Float64}(undef, nzval_nbr)
        new(
            dynamic_variables,
            dynamic_variables2,
            exogenous_variables,
            residuals,
            derivatives,
            temporary_values,
            nzval
        )
    end
end

function DynamicWs(context::Context)
    m = context.models[1]
    df = context.dynarefunctions
    dynamic_nbr = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr = sum(df.dynamic!.tmp_nbr[1:2])
    nzval_nbr = length(m.dynamic_g1_dynamic_rowval)
    return DynamicWs(m.endogenous_nbr, m.exogenous_nbr, dynamic_nbr, tmp_nbr, nzval_nbr)
end

struct StaticWs
    residuals::Vector{Float64}
    derivatives::Vector{Vector{Float64}}
    temporary_values::Vector{Float64}
    function StaticWs(endogenous_nbr::Int64, tmp_nbr::Int64)
        residuals = Vector{Float64}(undef, endogenous_nbr)
        derivatives = Vector{Vector{Float64}}(undef, 1)
        derivatives[1] = Vector{Float64}(undef, 0)
        temporary_values = Vector{Float64}(undef, tmp_nbr)
        new(residuals, derivatives, temporary_values)
    end
end

function StaticWs(context::Context)
    m = context.models[1]
    df = context.dynarefunctions
    tmp_nbr = sum(df.static!.tmp_nbr[1:2])
    return StaticWs(m.endogenous_nbr, tmp_nbr)
end

function inverse_order_of_dynare_decision_rule(m::Model)
    inverse_order_var = Vector{Int64}(undef, m.endogenous_nbr)
    for i = 1:m.n_static
        inverse_order_var[m.i_static[i]] = i
    end

    offset = m.n_static
    for i = 1:m.n_bkwrd
        inverse_order_var[m.i_bkwrd[i]] = i + offset
    end

    offset += m.n_bkwrd
    for i = 1:m.n_both
        inverse_order_var[m.i_both[i]] = i + offset
    end

    offset += m.n_both
    for i = 1:m.n_fwrd
        inverse_order_var[m.i_fwrd[i]] = i + offset
    end

    inverse_order_states = sortperm(cat(m.i_bkwrd, m.i_both; dims = 1))

    (inverse_order_var, inverse_order_states)
end

function get_initial_dynamic_endogenous_variables!(
    y::AbstractVector{Float64},
    data::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    lli::Matrix{Int64},
    period::Int64,
)
    m, n = size(lli)
    p = (period - 2) * n
    for j = 1:n
        k = lli[1, j]
        if k > 0
            y[k] = initialvalues[p+j]
        end
    end
    @inbounds for i = 2:m
        for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p+j]
            end
        end
        p += n
    end
end

function get_terminal_dynamic_endogenous_variables!(
    y::AbstractVector{Float64},
    data::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    lli::Matrix{Int64},
    period::Int64,
)
    m, n = size(lli)
    p = (period - 2) * n
    @inbounds for i = 1:m-1
        for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p+j]
            end
        end
        p += n
    end
    for j = 1:n
        k = lli[m, j]
        if k > 0
            y[k] = terminalvalues[j]
        end
    end
end

function get_dynamic_endogenous_variables!(
    y::AbstractVector{Float64},
    data::AbstractVector{Float64},
    lli::Matrix{Int64},
    period::Int64,
)
    m, n = size(lli)
    p = (period - 2) * n
    @inbounds for i = 1:m
        for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p+j]
            end
        end
        p += n
    end
end

function get_dynamic_endogenous_variables2!(
    y::AbstractVector{Float64},
    data::AbstractVector{Float64},
    endogenous_nbr,
    period::Int64,
)
    p = (period - 2) * endogenous_nbr
    y .= view(data, p .+ (1:3*endogenous_nbr))
end


"""
`get_dynamic_endogenous_variables!`(y::Vector{Float64}, data::Vector{Float64}, lli::Matrix{Int64})

sets the vector of dynamic variables ``y``, evaluated at the same values 
for all leads and lags and taken in ``data`` vector 
"""
function get_dynamic_endogenous_variables!(
    y::Vector{Float64},
    data::AbstractVector{Float64},
    lli::Matrix{Int64},
)
    for i = 1:size(lli, 2)
        value = data[i]
        for j = 1:size(lli, 1)
            k = lli[j, i]
            if k > 0
                y[k] = value
            end
        end
    end
end

function get_dynamic_endogenous_variables2!(
    y::Vector{Float64},
    data::AbstractVector{Float64},
    endogenous_nbr::Int64,
)
    @views begin
        k = 1:endogenous_nbr
        @show k
        @show data
        @show y[k]
        y[k] .= data
        y[endogenous_nbr .+ k] .= data
        y[2*endogenous_nbr .+ k] .= data
    end
    return y
end

"""
`get_dynamic_endogenous_variables!`(y::Vector{Float64}, data::Matrix{Float64}, lli::Matrix{Int64}, m::Model, period::Int64)

sets the vector of dynamic variables ``y`` with values in as many rows of ``data`` matrix
as there are leads and lags in the model. ``period`` is the current period.
"""
function get_dynamic_endogenous_variables!(
    y::Vector{Float64},
    data::AbstractMatrix{Float64},
    lli::Matrix{Int64},
    m::Model,
    period::Int64,
)
    for i = 1:size(lli, 2)
        p = period - m.maximum_lag - 1
        for j = 1:size(lli, 1)
            k = lli[j, i]
            if k > 0
                y[k] = data[p+j, i]
            end
        end
    end
end

"""
`get_dynamic_endogenous_variables!`(y::Vector{Float64}, data::Vector{Float64}, lli::Matrix{Int64}, m::Model, period::Int64)

sets the vector of dynamic variables ``y`` with values in as many rows of ``data`` matrix
as there are leads and lags in the model. ``period`` is the current period.
"""
function get_dynamic_endogenous_variables!(
    y::Vector{Float64},
    data::AbstractVector{Float64},
    lli::Matrix{Int64},
    m::Model,
    period::Int64,
)
    n = m.endogenous_nbr
    p = (period - m.maximum_lag - 1) * n
    for j = 1:size(lli, 1)
        for i = 1:n
            k = lli[j, i]
            if k > 0
                y[k] = data[p+i]
            end
        end
        p += n
    end
end

function get_exogenous_matrix(x::Vector{Float64}, exogenous_nbr::Int64)
    @debug "any(isnan.(x))=$(any(isnan.(x))) "
    x1 = reshape(x, Int(length(x) / exogenous_nbr), exogenous_nbr)
    @debug "any(isnan.(x1))=$(any(isnan.(x1))) "
    return x1
end

function get_dynamic_variables!(
    ws::DynamicWs,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    m::Model,
    period::Int64,
)
    lli = m.lead_lag_incidence
    @show endogenous
    get_dynamic_endogenous_variables!(ws.dynamic_variables, endogenous, lli)
    lx = length(ws.exogenous_variables)
    nrx = period + m.maximum_exo_lead
    required_lx = nrx * m.exogenous_nbr
    if lx < required_lx
        resize!(ws.exogenous_variables, required_lx)
        lx = required_lx
    end
    @debug "any(isnan.(ws.exognoues_variables))=$(any(isnan.(ws.exogenous_variables)))"
    x = get_exogenous_matrix(ws.exogenous_variables, m.exogenous_nbr)
    x .= transpose(exogenous)
end

function get_dynamic_residuals!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    x = get_dynamic_variables!(ws, endegenous, exogenous, m, period)
    Base.invokelatest(
        df.dynamic!.dynamic!,
        ws.temporary_values,
        ws.residuals,
        ws.dynamic_variables,
        x,
        params,
        steadystate,
        period,
    )
    return ws.residuals
end

function get_static_residuals!(
    ws::StaticWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    df::DynareFunctions,
)
    Base.invokelatest(
        df.static!.static!,
        ws.temporary_values,
        ws.residuals,
        endogenous,
        exogenous,
        params,
    )
    return ws.residuals
end

function dynamic_jacobian_matrix(ws::DynamicWs, m::Model)
    dynamic_nbr = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    vecjacobian =
        resize!(ws.derivatives[1], m.endogenous_nbr * (dynamic_nbr + m.exogenous_nbr))
    jacobian = reshape(vecjacobian, m.endogenous_nbr, dynamic_nbr + m.exogenous_nbr)
    fill!(jacobian, 0.0)
    return jacobian
end

function static_jacobian_matrix(ws::StaticWs, n::Int64)
    vecjacobian = resize!(ws.derivatives[1], n * n)
    jacobian = reshape(vecjacobian, n, n)
    fill!(jacobian, 0.0)
    return jacobian
end

"""
`get_dynamic_jacobian!`(ws::DynamicWs, params::Vector{Float64}, endogenous::AbstractVector{Float64}, exogenous::Vector{Float64}, m::Model, period::Int64)

sets the dynamic Jacobian matrix ``work.jacobian``, evaluated at ``endogenous`` and ``exogenous`` values, identical for all leads and lags
"""
function get_dynamic_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    x = get_dynamic_variables!(ws, endogenous, exogenous, m, period)
    jacobian = dynamic_jacobian_matrix(ws, m)
    Base.invokelatest(
        df.dynamic!.dynamic!,
        ws.temporary_values,
        ws.residuals,
        jacobian,
        ws.dynamic_variables,
        x,
        params,
        steadystate,
        period,
    )
    return jacobian
end

function get_initial_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    steadystate::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    lli = m.lead_lag_incidence
    get_initial_dynamic_endogenous_variables!(
        ws.dynamic_variables,
        endogenous,
        initialvalues,
        lli,
        period,
    )
    jacobian = dynamic_jacobian_matrix(ws, m)
    fill!(jacobian, 0.0)
    Base.invokelatest(
        df.dynamic!.dynamic!,
        ws.temporary_values,
        ws.residuals,
        jacobian,
        ws.dynamic_variables,
        exogenous,
        params,
        steadystate,
        period,
    )
    return jacobian
end

function get_terminal_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    steadystate::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    lli = m.lead_lag_incidence
    get_terminal_dynamic_endogenous_variables!(
        ws.dynamic_variables,
        endogenous,
        terminalvalues,
        lli,
        period,
    )
    jacobian = dynamic_jacobian_matrix(ws, m)
    Base.invokelatest(
        df.dynamic!.dynamic!,
        ws.temporary_values,
        ws.residuals,
        jacobian,
        ws.dynamic_variables,
        exogenous,
        params,
        steadystate,
        period,
    )
    return jacobian
end


"""
`get_dynamic_jacobian!`(ws::DynamicWs,
                        params::Vector{Float64},
                        endogenous::AbstractVecOrMat{Float64},
                        exogenous::Matrix{Float64},
                        steadystate::Vector{Float64},
                        m::Model,
                        df::DynareFunctions,
                        period::Int64
                       )

sets the dynamic Jacobian matrix ``ws.jacobian``, evaluated with ``endogenous`` and ``exogenous`` values taken
around ``period`` 
"""
function get_dynamic_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVecOrMat{Float64},
    exogenous::Matrix{Float64},
    steadystate::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    lli = m.lead_lag_incidence
    @show length(endogenous)
    get_dynamic_endogenous_variables!(ws.dynamic_variables, endogenous, lli, m, period)
    dynamic_nbr = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    jacobian = dynamic_jacobian_matrix(ws, m)
    y = ws.dynamic_variables
    x = exogenous
    it_ = period
    Base.invokelatest(
        df.dynamic!.dynamic!,
        ws.temporary_values,
        ws.residuals,
        jacobian,
        ws.dynamic_variables,
        exogenous,
        params,
        steadystate,
        period,
    )
    return jacobian
end

"""
`get_static_jacobian!`(ws::StaticWs, params::Vector{Float64}, endogenous::AbstractVector{Float64}, exogenous::Vector{Float64})

sets the static Jacobian matrix ``work.jacobian``, evaluated at ``endogenous`` and ``exogenous`` values
"""
function get_static_jacobian!(
    ws::StaticWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    m::Model,
    df::DynareFunctions,
)
    @debug "any(isnan.(exognous))=$(any(isnan.(exogenous)))"
    jacobian = static_jacobian_matrix(ws, m.endogenous_nbr)
    Base.invokelatest(
        df.static!.static!,
        ws.temporary_values,
        ws.residuals,
        jacobian,
        endogenous,
        exogenous,
        params,
    )
    return jacobian
end

function get_abc!(
    a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64},
    jacobian::AbstractMatrix{Float64},
    m::Model,
)
    i_rows = (m.n_static+1):m.endogenous_nbr
    fill!(a, 0.0)
    fill!(b, 0.0)
    fill!(c, 0.0)
    ws.a[:, ws.forward_indices_d] .=
        view(jacobian, i_rows, ws.backward_nbr .+ ws.current_nbr .+ (1:ws.forward_nbr))
    ws.b[:, ws.current_dynamic_indices_d] .=
        view(jacobian, i_rows, ws.backward_nbr .+ ws.current_dynamic_indices)
    ws.c[:, ws.backward_indices_d] .= view(jacobian, i_rows, 1:ws.backward_nbr)
end

function get_de!(ws::LinearRationalExpectationsWs, jacobian::AbstractMatrix{Float64})
    n1 = ws.backward_nbr + ws.forward_nbr - ws.both_nbr
    fill!(ws.d, 0.0)
    fill!(ws.e, 0.0)
    i_rows = (ws.static_nbr+1):ws.endogenous_nbr
    ws.d[1:n1, ws.icolsD] .= jacobian[i_rows, ws.jcolsD]
    ws.e[1:n1, ws.icolsE] .= -jacobian[i_rows, ws.jcolsE]
    u = Matrix{Float64}(I, ws.both_nbr, ws.both_nbr)
    i_rows = n1 .+ (1:ws.both_nbr)
    ws.d[i_rows, ws.colsUD] .= u
    ws.e[i_rows, ws.colsUE] .= u
end

function get_sparse_dynamic_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVecOrMat{Float64},
    exogenous::Matrix{Float64},
    steady_state::Vector{Float64},
    m::Model,
    df::DynareFunctions,
    period::Int64,
)
    get_dynamic_endogenous_variables2!(ws.dynamic_variables2, endogenous, m.endogenous_nbr)
    global SparseDynamicResidTT! = df.SparseDynamicResidTT!
    global SparseDynamicG1TT! = df.SparseDynamicG1TT!
    df.SparseDynamicG1!(ws.temporary_values,
                        ws.nzval,
                        ws.dynamic_variables2,
                        exogenous,
                        params,
                        steady_state,
                        length(ws.temporary_values) > 0)
    @show m.endogenous_nbr
    @show 3*m.endogenous_nbr + m.exogenous_nbr
    @show m.dynamic_g1_sparse_colptr
    @show m.dynamic_g1_sparse_rowval
    @show ws.nzval
    jacobian = SparseMatrixCSC{Float64, Int64}(m.endogenous_nbr,
                                               3*m.endogenous_nbr + m.exogenous_nbr,
                                               m.dynamic_g1_sparse_colptr,
                                               m.dynamic_g1_sparse_rowval,
                                               ws.nzval)
    return jacobian
end

nearbyint(x::T) where T <: Real  = (abs((x)-floor(x)) < abs((x)-ceil(x)) ? floor(x) : ceil(x))

function get_power_deriv(x::T, p::T, k::Int64) where T <: Real
    if (abs(x) < 1e-12 && p > 0 && k > p && abs(p-nearbyint(p)) < 1e-12 )
        return 0.0
    else
        dxp = x^(p-k)
        for i = 1:k
	     dxp *= p
	     p -= 1
	 end
	 return dxp
    end
end

