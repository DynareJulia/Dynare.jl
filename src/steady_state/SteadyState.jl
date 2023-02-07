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
        homotopy_mode = None
        homotopy_steps = 10
        nocheck = false
        for (k, v) in pairs(options)
            if k == "noprint"
                display = false
            elseif k == "maxit" && v::Bool
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
            elseif k == "steadystate.nocheck" && v::Bool
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
    steady!(context::Context, field::Dict{String, Any})

Compute the steady state of the model and set the result in `context`
"""
function steady!(context::Context, field::Dict{String,Any})
    model = context.models[1]
    modfileinfo = context.modfileinfo
    options = SteadyOptions(get(field, "options", field))
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
    compute_steady_state!(context, field)
    if options.display
        if isempty(trends.endogenous_terminal_steady_state)
            steadystate_display(trends.endogenous_steady_state, context)
        else
            steadystate_display(trends.endogenous_steady_state, context, title_complement = "Initial")
            steadystate_display(trends.endogenous_terminal_steady_state, context, title_complement = "Terminal")
        end
    end
end

function compute_steady_state!(context::Context, field::Dict{String,Any})
    options = SteadyOptions(get(field, "options", field))
    model = context.models[1]
    modfileinfo = context.modfileinfo
    trends = context.results.model_results[1].trends
    work = context.work
    endogenous_nbr = context.models[1].endogenous_nbr
    if (
        modfileinfo.has_steadystate_file &&
        length(work.analytical_steadystate_variables) ==
        endogenous_nbr
    )
        evaluate_steady_state!(trends.endogenous_steady_state,
                               trends.exogenous_steady_state,
                               work.params)
        !options.nocheck && check_steady_state(StaticWs(context),
                                               trends.endogenous_steady_state,
                                               trends.exogenous_steady_state,
                                               work.params,
                                               options.tolf)
        if !isempty(trends.exogenous_terminal_steady_state)
            evaluate_steady_state!(trends.endogenous_terminal_steady_state,
                                   trends.exogenous_terminal_steady_state,
                                   work.params)
            check_steady_state(StaticWs(context),
                               trends.endogenous_terminal_steady_state,
                               trends.exogenous_terminal_steady_state,
                               work.params,
                               options.tolf)
        end
    elseif context.modfileinfo.has_ramsey_model
        x0 = zeros(model.endogenous_nbr)
        !isempty(trends.endogenous_steady_state) &&
            (x0 .= Float64.(trends.endogenous_steady_state))
        solve_ramsey_steady_state!(context, x0, options)
    else
        # initial steady state
        exogenous = zeros(model.exogenous_nbr)
        !isempty(trends.exogenous_steady_state) &&
            (exogenous .= Float64.(trends.exogenous_steady_state))
        x0 = zeros(model.endogenous_nbr)
        !isempty(trends.endogenous_steady_state) &&
            (x0 .= Float64.(trends.endogenous_steady_state))
        trends.endogenous_steady_state .=
            solve_steady_state!(context, x0, exogenous, options)
        # terminal steady state
        if modfileinfo.has_endval
            !isempty(trends.exogenous_terminal_steady_state) &&
                (exogenous .= Float64.(trends.exogenous_terminal_steady_state))
            x0 = zeros(model.endogenous_nbr)
            !isempty(trends.endogenous_terminal_steady_state) &&
                (x0 .= Float64.(trends.endogenous_terminal_steady_state))
            trends.endogenous_terminal_steady_state .=
                solve_steady_state!(context, x0, exogenous, options)
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
                        exogenous::Vector{<:Real},
                        options::SteadyOptions)

Solve the steady state numerically
"""
function solve_steady_state!(context::Context,
                             x0::AbstractVector{<:Real},
                             exogenous::AbstractVector{<:Real},
                             options::SteadyOptions)
        try
            return solve_steady_state_!(context, x0, exogenous, options)
        catch e
            if !isempty(x0) && isa(e, Dynare.DynareSteadyStateComputationFailed)
                i = 1
                while i <= 5
                    x00 = rand(0.95:0.01:1.05, length(x0)).*x0
                    try
                        return solve_steady_state_!(context, x00, exogenous, options)
                    catch
                    end
                    i += 1
                end
                if i > options.maxit
                    rethrow(e)
                end
            else
                rethrow(e)
            end
        end        
end

"""
    solve_steady_state_!(context::Context,
                         x0::Vector{Float64},
                         exogenous::AbstractVector{<:Real})

Solve the static model to obtain the steady state
"""
function solve_steady_state_!(context::Context,
                              x0::AbstractVector{<:Real},
                              exogenous::AbstractVector{<:Real},
                              options)
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
    
    J! = make_static_jacobian(ws.temporary_values,
                              exogenous,
                              work.params)

    results = context.results.model_results[1]

    function f!(residuals::AbstractVector{Float64}, x::AbstractVector{Float64})
        DFunctions.static!(ws.temporary_values, residuals, x, exogenous, work.params)
    end

    of = OnceDifferentiable(f!, J!, vec(x0), residuals, A)
    result = nlsolve(of, x0; method = :robust_trust_region, show_trace = true, ftol = options.tolf)
    @debug result
    if converged(result)
        return result.zero
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed())
    end
end

"""
    solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}, options)

Solve numerically for the steady state of a Ramsey problem
"""
function solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}, options)
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
    
    f!(x00) <  options.tolf && copy!(endogenous_steady_state, endogenous)
    if unknown_variable_nbr == 0
        @debug "Steady state computation failed"
        throw(DynareSteadyStateComputationFailed())
    else
        result = optimize(f!, x00, LBFGS(), Optim.Options(f_tol=1e-6))
        @debug result
        if Optim.converged(result) && abs(Optim.minimum(result)) < options.tolf
            view(endogenous, unknown_variable_indices) .= Optim.minimizer(result)
            context.modfileinfo.has_auxiliary_variables &&
                DFunctions.static_auxiliary_variables!(endogenous, exogenous, params)
            context.modfileinfo.has_steadystate_file &&
                DFunctions.steady_state!(endogenous, exogenous, params)
            # get Lagrance multipliers
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
    f!(residuals, x) = DFunctions.static!(temp_val,
                                          residuals,
                                          x,
                                          exogenous,
                                          params)
    return f!
end

function make_static_jacobian(temp_val::AbstractVector{T},
                              exogenous::AbstractVector{T},
                              params::AbstractVector{T}) where T <: Real
    f!(A, x) = DFunctions.static_derivatives!(temp_val,
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

    
