using Optim

struct DynareSteadyStateComputationFailed <: Exception end
Base.showerror(io::IO, e::DynareSteadyStateComputationFailed) = print("""
                                                                      Dynare couldn't compute the steady state.
                                                                      Either there is no solution or the guess values
                                                                      are too far from the solution
                                                                       """)

"""
NonLinearSolveAlgos - enumerate

Nonlinear system of equations available algorithms:
- trustregion: trust-region algorithm from NLsolve package
"""
@enum NonLinearSolveAlgos trustregion

"""
HomotopyModes - enumerate
Homotopy modes:
- None
- SimultaneousFixedSteps: move all declared parameters in a fixed number of steps
- SingleParameterFixedSteps: move a single parameter in a fixed number of steps
- SimultaneousAdaptative: move all declared parameters in a adaptative manner
"""
@enum HomotopyModes None SimultaneousFixedSteps SingleParameterFixedSteps SimultaneousAdaptative

"""
SteadyOptions type
    display::Bool - whether to display results [true]
    maxit::Int64 - maximum number of iterations [50]
    tolf::Float64 - tolerance criterium for equations error [1e-8]
    tolx::Float64 - norm difference in x between two successive iterates under which convergence is declared.  [0]
    solve_algo::NonLinearSolveAlgos - algorithm for nonlinear equations solver [trustregion]
    homotopy_mode::HomotopyModes - homotopy mode [None] 
    homotopy_steps::Int64 - homotopy steps [10]
    nocheck::Bool - don't check steady state values provided by the user [false]
"""
struct SteadyOptions
    display::Bool
    maxit::Int64
    tolf::Float64
    tolx::Float64
    solve_algo::NonLinearSolveAlgos
    homotopy_mode::HomotopyModes
    homotopy_steps::Int64
    nocheck::Bool
    function SteadyOptions(options::Dict{String,Any})
        display = true
        maxit = 50
        tolf = cbrt(eps())
        tolx = 0.0
        solve_algo = trustregion
        homotopy_mode = SimultaneousFixedSteps
        homotopy_steps = 0
        nocheck = false
        for (k, v) in pairs(options)
            if k == "noprint"
                display = false
            elseif k == "steady.maxit"
                maxit = v::Int64
            elseif k == "solve_tolf"
                tolf = v::Float64
            elseif k == "solve_tolx"
                tolx = v::Float64
            elseif k == "solve_algo"
                solve_algo = v::NonLinearSolveAlgos
            elseif k == "homotopy_mode"
                homotopy_mode = v::HomotopyModes
            elseif k == "homotopy_steps"
                homotopy_steps = v::Int64
            elseif k == "steadystate.nocheck"
                nocheck = true
            end
        end
        new(display, maxit, tolf, tolx, solve_algo, homotopy_mode, homotopy_steps, nocheck)
    end
end

function set_or_zero!(x, a, n)
    resize!(x, n)
    if !isempty(a)
        copyto!(x,a)
    else
        fill!(x, 0.0)
    end
end
        
"""
    steadystate!(; context::Context=context, display::Bool = true,
            maxit::Int = 50, nocheck::Bool = false, tolf::Float64 = cbrt(eps()),
            tolx::Float64 = 0.0)
                             
 # Keyword arguments
 - `context::Context=context`: context in which the simulation is computed
 - `homotopy_steps::Int=0`: number of homotopy steps
 - `display::Bool=true`: whether to display the results
 - `maxit::Int=50` maximum number of iterations
 - `nocheck::Bool=false`: don't check the steady state
 - `tolf::Float64=cbrt(eps())`: tolerance for the norm of residualts
 - `tolx::Float64=0`: tolerance for the norm of the change in the result
"""
function steadystate!(; context::Context=context,
                      display = true,
                      homotopy_mode = SimultaneousFixedSteps,
                      homotopy_steps = 0,
                      maxit = 50,
                      nocheck = false,
                      nonlinear_solve_algo = TrustRegion(),
                      tolf = cbrt(eps()),
                      tolx = 0.0
                      )
    model = context.models[1]
    modfileinfo = context.modfileinfo
    trends = context.results.model_results[1].trends
    work = context.work
    
    set_or_zero!(trends.endogenous_steady_state,
                 work.initval_endogenous,
                 model.endogenous_nbr)
    set_or_zero!(trends.exogenous_steady_state,
                 work.initval_exogenous,
                 model.exogenous_nbr)
    if modfileinfo.has_endval
        set_or_zero!(trends.endogenous_terminal_steady_state,
                     work.endval_endogenous,
                     model.endogenous_nbr)
        set_or_zero!(trends.exogenous_terminal_steady_state,
                     work.endval_exogenous,
                     model.exogenous_nbr)
    end
    if homotopy_steps == 0
        compute_steady_state!(context, maxit = maxit, nocheck = nocheck, nonlinear_solve_algo = nonlinear_solve_algo, tolf = tolf)
    else
        homotopy_steady!(context, homotopy_mode, homotopy_steps, maxit = maxit, nonlinear_solve_algo, tolf = tolf)
    end
    if display
        if isempty(trends.endogenous_terminal_steady_state)
            steadystate_display(trends.endogenous_steady_state, context)
        else
            steadystate_display(trends.endogenous_steady_state, context, title_complement = "Initial")
            steadystate_display(trends.endogenous_terminal_steady_state, context, title_complement = "Terminal")
        end
    end
end

function steady_!(context, field::Dict{String, Any})
    options = SteadyOptions(get(field, "options", field))
    steadystate!(context = context,
                 display = options.display,
                 homotopy_mode = options.homotopy_mode,
                 homotopy_steps = options.homotopy_steps,
                 maxit = options.maxit,
                 nocheck = options.nocheck,
                 tolf = options.tolf,
                 tolx = options.tolx
                 )
    return nothing
end

function check_steadystate_model(context::Context, nocheck)
    if context.modfileinfo.has_steadystate_file
        if (length(context.work.analytical_steadystate_variables) ==
            context.models[1].endogenous_nbr)
            return true
        else
            if nocheck
                oem = context.models[1].original_endogenous_nbr
                endogenous_names = get_endogenous(context.symboltable)
                missing_indices = setdiff(1:oem, context.work.analytical_steadystate_variables)
                missing_variables = endogenous_names[missing_indices]
                error("Option nocheck is used but variables are missing in steady_state_model block: "*
                       join(missing_variables, ", "))
            end
            return false
        end
    else
        return false
    end 
end 

function compute_steady_state!(context::Context; maxit = 50, nocheck = false, nonlinear_solve_algo = TrustRegion(), tolf = cbrt(eps()))
    model = context.models[1]
    modfileinfo = context.modfileinfo
    trends = context.results.model_results[1].trends
    work = context.work
    endogenous_nbr = context.models[1].endogenous_nbr
    exogenous_nbr = context.models[1].exogenous_nbr
    isempty(trends.endogenous_steady_state) && (trends.endogenous_steady_state = Vector{Float64}(undef, endogenous_nbr))
    isempty(trends.exogenous_steady_state) && (trends.exogenous_steady_state = Vector{Float64}(undef, exogenous_nbr))
    if check_steadystate_model(context, nocheck)
        evaluate_steady_state!(trends.endogenous_steady_state,
                               trends.exogenous_steady_state,
                               work.params)
        !nocheck && check_steady_state(StaticWs(context),
                                       trends.endogenous_steady_state,
                                       trends.exogenous_steady_state,
                                       work.params, tolf)
        if !isempty(trends.exogenous_terminal_steady_state)
            evaluate_steady_state!(trends.endogenous_terminal_steady_state,
                                   trends.exogenous_terminal_steady_state,
                                   work.params)
            !nocheck && check_steady_state(StaticWs(context),
                                           trends.endogenous_terminal_steady_state,
                                           trends.exogenous_terminal_steady_state,
                                           work.params, tolf)
        end
    elseif context.modfileinfo.has_ramsey_model
        x0 = zeros(model.endogenous_nbr)
        !isempty(trends.endogenous_steady_state) &&
            (x0 .= Float64.(trends.endogenous_steady_state))
        solve_ramsey_steady_state!(context, x0, tolf = tolf)
    else
        # initial steady state
        exogenous = zeros(model.exogenous_nbr)
        !isempty(trends.exogenous_steady_state) &&
            (exogenous .= Float64.(trends.exogenous_steady_state))
        x0 = zeros(model.endogenous_nbr)
        !isempty(trends.endogenous_steady_state) &&
            (x0 .= Float64.(trends.endogenous_steady_state))
        if !nocheck
            trends.endogenous_steady_state .=
                solve_steady_state!(context, x0, exogenous, maxit = maxit, nonlinear_solve_algo = nonlinear_solve_algo, tolf = tolf)
        end
        # terminal steady state
        if modfileinfo.has_endval
            !isempty(trends.exogenous_terminal_steady_state) &&
                (exogenous .= Float64.(trends.exogenous_terminal_steady_state))
            x0 = zeros(model.endogenous_nbr)
            !isempty(trends.endogenous_terminal_steady_state) &&
                (x0 .= Float64.(trends.endogenous_terminal_steady_state))
            !nocheck && (trends.endogenous_terminal_steady_state .=
                solve_steady_state!(context, x0, exogenous, maxit = maxit, nonlinear_solve_algo = nonlinear_solve_algo, tolf = tolf))
        end
    end
end

    
"""
    steadystate_display(steady_state, context; title_complement = "")

Display the steady state of the model
"""
function steadystate_display(steady_state::AbstractVector{<:Real},
                             context::Context;
                             title_complement = "")
    m = context.models[1]
    results = context.results.model_results[1]
    endogenous_names = get_endogenous_longname(context.symboltable)
    n = m.original_endogenous_nbr
    labels = endogenous_names[1:n]
    data = Matrix{Any}(undef, n, 2)
    for i = 1:n
        data[i, 1] = labels[i]
        data[i, 2] = steady_state[i]
    end
    if isempty(title_complement)
        title = "Steady state"
    else
        title = "$title_complement steady state"
    end
    println(" ")
    dynare_table(data, title, columnheader = false)
end

"""
    evaluate_steady_state!(endogenous::Vector{<: Real},
                           exogenous::Vector{<: Real},
                           params::AbstractVector{<: Real})

Evaluate the steady state function provided by the user.
"""
function evaluate_steady_state!(
    endogenous::Vector{<: Real},
    exogenous::Vector{<: Real},
    params::AbstractVector{<: Real},
)
    DFunctions.steady_state!(
        endogenous,
        exogenous,
        params,
    )
end

"""
    solve_steady_state!(context::Context,
                        x0::Vector{<:Real},
                        exogenous::Vector{<:Real};
                        maxit = 50,
                        tolf = cbrt(eps()))

Solve the steady state numerically
"""
function solve_steady_state!(context::Context,
                             x0::AbstractVector{<:Real},
                             exogenous::AbstractVector{<:Real};
                             maxit = 50,
                             nonlinear_solve_algo = TrustRegion(),
                             tolf = cbrt(eps()))
    try
            return solve_steady_state_!(context, x0, exogenous; maxit = maxit, nonlinear_solve_algo = nonlinear_solve_algo, tolf = tolf)
        catch e
            if !isempty(x0) && isa(e, Dynare.DynareSteadyStateComputationFailed)
                i = 1
                while i <= 0
                    x00 = rand(0.95:0.01:1.05, length(x0)).*x0
                    try
                        return solve_steady_state_!(context, x00, exogenous, maxit = maxit, nonlinear_solve_algo = nonlinear_solve_algo, tolf = tolf)
                    catch
                    end
                    i += 1
                end
                if i > maxit
                    rethrow(e)
                end
            else
                rethrow(e)
            end
        end        
end
using LinearSolve
"""
    solve_steady_state_!(context::Context;
                         x0::Vector{Float64},
                         exogenous::AbstractVector{<:Real};
                         maxit = 50,
                         nonlinear_solve_algo = TrustRegion(),
                         show_trace = false,
                         tolf = cbrt(eps()))

Solve the static model to obtain the steady state
"""
function solve_steady_state_!(context::Context,
                              x0::AbstractVector{<:Real},
                              exogenous::AbstractVector{<:Real};
                              maxit = 50,
                              nonlinear_solve_algo = TrustRegion(),
                              show_trace = false,
                              tolf = cbrt(eps()))
    ws = StaticWs(context)
    model = context.models[1]
    work = context.work
    results = context.results.model_results[1]
    params = work.params
    A = ws.derivatives[1]
    residuals = zeros(model.endogenous_nbr)

    f! = make_static_residuals(ws.temporary_values,
                               exogenous,
                               work.params)
    
    j! = make_static_jacobian(ws.temporary_values,
                              exogenous,
                              work.params)

    results = context.results.model_results[1]
    j!(A, x0, params)
    lp = LinearProblem(A, residuals)
    dy1 = LinearSolve.solve(lp)
    dy2 = LinearSolve.solve(lp, PardisoJL())

    # of = OnceDifferentiable(f!, J!, vec(x0), residuals, A)
    # result = nlsolve(of, x0; method = :robust_trust_region, show_trace = false, ftol = tolf, iterations = maxit)
    fj = NonlinearFunction(f!, jac = j!, jac_prototype = A)
    prob = NonlinearProblem(fj, x0, params)
    
    result = NonlinearSolve.solve(prob, nonlinear_solve_algo, show_trace=Val(show_trace), abstol = tolf)
    @debug result
    if result.retcode == ReturnCode.Success
        return result.u
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed())
    end
end

"""
    solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}; tolf = cbrt(eps()) )

Solve numerically for the steady state of a Ramsey problem
"""
function solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}; tolf = cbrt(eps()))
    ws = StaticWs(context)
    model = context.models[1]
    work = context.work
    params = work.params
    endogenous = zeros(model.endogenous_nbr)
    results = context.results.model_results[1]
    endogenous_steady_state = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    orig_endo_nbr = model.original_endogenous_nbr
    mult_indices =
        orig_endo_nbr .+ findall(a -> a["type"] == 6, context.models[1].auxiliary_variables)
    mult_nbr = length(mult_indices)
    mult = zeros(mult_nbr)
    orig_endo_aux_nbr = mult_indices[1] - 1
    unknown_variable_indices =
        setdiff!(collect(1:model.original_endogenous_nbr), work.analytical_steadystate_variables)
    unknown_variable_nbr = length(unknown_variable_indices)
    M = zeros(orig_endo_nbr, mult_nbr)
    U1 = zeros(orig_endo_nbr)
    residuals = zeros(model.endogenous_nbr)
    x00 = zeros(unknown_variable_nbr)
    x00 .= view(x0, unknown_variable_indices)

    f! = make_static_ramsey_residuals(ws.temporary_values,
                                      endogenous,
                                      exogenous,
                                      work.params,
                                      U1,
                                      M,
                                      mult_indices,
                                      unknown_variable_indices,
                                      orig_endo_nbr,
                                      ws)
    
    if f!(x00) <  tolf
        copy!(results.trends.endogenous_steady_state, endogenous)
        return
    end 
    if unknown_variable_nbr == 0 
        @debug "Steady state computation failed"
        throw(DynareSteadyStateComputationFailed())
    else
        result = Optim.optimize(f!, x00, LBFGS(), Optim.Options(f_tol=1e-6))
        if Optim.converged(result) && abs(Optim.minimum(result)) < tolf
            view(endogenous, unknown_variable_indices) .= Optim.minimizer(result)
            context.modfileinfo.has_auxiliary_variables &&
                DFunctions.static_auxiliary_variables!(endogenous, exogenous, params)
            context.modfileinfo.has_steadystate_file &&
                DFunctions.steady_state!(endogenous, exogenous, params)
            # get Lagrange multipliers
            f!(Optim.minimizer(result))
            results.trends.endogenous_steady_state .= endogenous
        else
            @debug "Steady state computation failed with\n $result"
            throw(DynareSteadyStateComputationFailed())
        end
    end
end

function check_steady_state(ws, endogenous, exogenous, params, tolf)
    residuals = get_static_residuals!(ws, params, endogenous, exogenous)
    if norm(residuals) > tolf
        display_residuals(residuals)
        throw(ErrorException("The steady_state_model block doesn't solve the static model!"))
    end
end

function display_residuals(residuals)
    println("\nThe steady_state_model block does'nt solve the static model")
    println("Equation Residuals")
    for (i, r) in enumerate(residuals)
        abs(r) > eps() && println("$i       $r")
    end
end

function make_static_residuals(temp_val::AbstractVector{T},
                               exogenous::AbstractVector{T},
                               params::AbstractVector{T}) where T <: Real
    function f!(residuals, x, params)
        try
            DFunctions.static!(temp_val,
                               residuals,
                               x,
                               exogenous,
                               params)
        catch
            return Inf*ones(length(residuals))
        end
    end
    
    return f!
end

function make_static_jacobian(temp_val::AbstractVector{T},
                              exogenous::AbstractVector{T},
                              params::AbstractVector{T}) where T <: Real
    f!(A, x, params) = DFunctions.static_derivatives!(temp_val,
                                              A,
                                              x,
                                              exogenous,
                                              params)
    return f!
end
    
function make_static_ramsey_residuals(temp_val::AbstractVector{T},
                                      endogenous::AbstractVector{T},
                                      exogenous::AbstractVector{T},
                                      params::AbstractVector{T},
                                      U1::AbstractVector{T},
                                      M::AbstractMatrix{T},
                                      mult_indices::AbstractVector{N},
                                      unknown_variable_indices::AbstractVector{N},
                                      orig_endo_nbr::N,
                                      ws::StaticWs
                                      ) where {T <: Real, N <: Integer}

    function f!(x::AbstractVector{Float64})
        view(endogenous, unknown_variable_indices) .= x
        view(endogenous, mult_indices) .= 0
        # Lagrange multipliers are kept to zero
        context.modfileinfo.has_auxiliary_variables &&
            DFunctions.static_auxiliary_variables!(endogenous, exogenous, params)
        context.modfileinfo.has_steadystate_file &&
            DFunctions.steady_state!(endogenous, exogenous, params)
        residuals = get_static_residuals!(ws, params, endogenous, exogenous)
        U1 .= view(residuals, 1:orig_endo_nbr)
        A = get_static_jacobian!(ws, params, endogenous, exogenous, context.models[1])
        M .= view(A, 1:orig_endo_nbr, mult_indices)
        mult = view(endogenous, mult_indices)
        mult .= -M \ U1
        view(residuals, 1:orig_endo_nbr) .= U1 .+ M * mult
        res1 = sum(x-> x*x, residuals)
        res2 = norm(residuals)

        return res2
    end

    return f!
end

function make_partial_static_residuals(temp_val::AbstractVector{T},
    exogenous::AbstractVector{T},
    params::AbstractVector{T},
    unknown_variable_indices::AbstractVector{N}) where {T <: Real, N <: Integer}

    function f!(residuals, x)
        view(endogenous, unknown_variable_indices) .= x
        context.modfileinfo.has_auxiliary_variables &&
            DFunctions.static_auxiliary_variables!(endogenous, exogenous, params)
        context.modfileinfo.has_steadystate_file &&
            DFunctions.steady_state!(endogenous, exogenous, params)
        DFunctions.static!(temp_val,
            residuals,
            endogenous,
            exogenous,
            params)
        return norm(residuals)
    end

    return f!
end

function homotopy_steady!(context::Context, homotopy_mode::HomotopyModes, homotopy_steps::Int64; maxit = 50, tolf = cbrt(eps()))
    hs = context.work.homotopy_setup
    initval_endogenous = context.work.initval_endogenous
    initval_exogenous = context.work.initval_exogenous
    params = context.work.params
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    exogenous_steadystate = context.results.model_results[1].trends.exogenous_steady_state
    if homotopy_mode == SimultaneousFixedSteps
        length(steadystate) > 0 && (initval_endogenous .= steadystate) 
        steps = zeros(length(hs))
        for (i,v) in enumerate(hs)
            if v.type == Dynare.Parameter
                startvalue = (ismissing(v.startvalue)) ? params[v.index] : v.startvalue
                steps[i] = Float64((v.endvalue - startvalue)/homotopy_steps)
                params[v.index] = startvalue
            elseif v.type == Dynare.Exogenous
                startvalue = (ismissing(v.startvalue)) ? exogenous_steadystate[v.index] : v.startvalue
                steps[i] = (v.endvalue - startvalue)/homotopy_steps
                exogenous_steadystate[v.index] = startvalue
            end
        end
        compute_steady_state!(context, maxit = maxit, nocheck = false, tolf = tolf)
        for i = 1:homotopy_steps
            for (i, v) in enumerate(hs)
                v.type == Dynare.Exogenous && (exogenous_steadystate[v.index] += steps[i])
                v.type == Dynare.Parameter && (params[v.index] += steps[i])
            end
            # if steadystate exists, use it as guess values
            length(steadystate) > 0 && (initval_endogenous .= steadystate) 
            compute_steady_state!(context, maxit = maxit, nocheck = false, tolf = tolf)
        end
    else
        error("homotopy_mode = $homotopy_mode isn't implemented")
    end
    length(steadystate) > 0 && (initval_endogenous .= steadystate) 
    length(exogenous_steadystate) > 0 && (initval_exogenous .= exogenous_steadystate) 
end
