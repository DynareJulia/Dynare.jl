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
        tolf = 1e-5 #cbrt(eps())
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
            elseif k == "nocheck" && v::Bool
                nochecl = true
            end
        end
        new(display, maxit, tolf, tolx, solve_algo, homotopy_mode, homotopy_steps, nocheck)
    end
end

"""
    steady!(context::Context, field::Dict{String, Any})

Compute the steady state of the model and set the result in `context`
"""
function steady!(context::Context, field::Dict{String,Any})
    options = SteadyOptions(get(field, "options", field))
    compute_steady_state!(context, field)
    if options.display
        steadystate_display(context)
    end
end

function compute_steady_state!(context::Context, field::Dict{String,Any})
    options = SteadyOptions(get(field, "options", field))
    modfileinfo = context.modfileinfo
    results = context.results.model_results[1]
    work = context.work
    endogenous_nbr = context.models[1].endogenous_nbr
    if (
        modfileinfo.has_steadystate_file &&
        length(context.work.analytical_steadystate_variables) ==
        endogenous_nbr
    )
        evaluate_steady_state!(results, work.params)
    else
        initval_endogenous = work.initval_endogenous
        initval_exogenous = work.initval_exogenous
        # will fail if missing values are encountered
        if size(initval_endogenous, 1) > 0
            x0 = Float64.(vec(view(context.work.initval_endogenous, 1, :)))
        else
            x0 = zeros(endogenous_nbr)
        end
        if size(initval_exogenous, 1) > 0
            copy!(
                results.trends.exogenous_steady_state,
                Float64.(vec(view(context.work.initval_exogenous, 1, :))),
            )
        else
            fill!(results.trends.exogenous_steady_state, 0.0)
        end
        try
            solve_steady_state!(context, x0, options)
        catch e
            if length(x0) > 0 && isa(e, Dynare.DynareSteadyStateComputationFailed)
                i = 1
                while i <= 5
                    x00 = rand(0.95:0.01:1.05, length(x0)).*x0
                    try
                        solve_steady_state!(context, x00, options)
                        break
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
end

    
"""
    steadystate_display(context::Context)

Display the steady state of the model
"""
function steadystate_display(context::Context)
    m = context.models[1]
    results = context.results.model_results[1]
    endogenous_names = get_endogenous_longname(context.symboltable)
    n = m.original_endogenous_nbr
    steady_state = results.trends.endogenous_steady_state[1:n]
    labels = endogenous_names[1:n]
    data = Matrix{Any}(undef, n, 2)
    for i = 1:n
        data[i, 1] = labels[i]
        data[i, 2] = steady_state[i]
    end
    title = "Steady state"
    dynare_table(data, title, columnheader = false)
end

"""
    evaluate_steady_state!(results::ModelResults,
                           steady_state!::Function,
                           params::AbstractVector{Float64})

Evaluate the steady state function provided by the user.

The steady state is stored in `context.results.trends.endogenous_steady_state`
"""
function evaluate_steady_state!(
    results::ModelResults,
    params::AbstractVector{Float64},
)
    fill!(results.trends.exogenous_steady_state, 0.0)
    DFunctions.steady_state!(
        results.trends.endogenous_steady_state,
        results.trends.exogenous_steady_state,
        params,
    )
end

"""
    sparse_static_jacobian(ws, params, x, exogenous, m, df)

Return the sparse Jacobina of a static model
"""
function sparse_static_jacobian(ws, params, x, exogenous, m, df)
    J = get_static_jacobian!(ws, params, x, exogenous, m, df)
    return sparse(J)
end

"""
    solve_steady_state!(context::Context,
                        x0::Vector{Float64},
                        options::SteadyOptions)

Solve the steady state numerically
"""
function solve_steady_state!(context::Context,
                             x0::AbstractVector{Float64},
                             options::SteadyOptions)
    if context.modfileinfo.has_ramsey_model
        solve_ramsey_steady_state!(context, x0, options)
    else
        solve_steady_state_!(context, x0, options)
    end
end

"""
    solve_steady_state_!(context::Context,
                                 x0::Vector{Float64})

Solve the static model to obtain the steady state
"""
function solve_steady_state_!(context::Context, x0::AbstractVector{Float64}, options)
    ws = StaticWs(context)
    m = context.models[1]
    w = context.work
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state

    function J1!(A, x::AbstractVector{Float64})
        DFunctions.static_derivatives!(ws.temporary_values,
                                       A,
                                       x,
                                       exogenous,
                                       w.params)
    end
    A = ws.derivatives[1]
    solve_steady_state_core!(context, x0, J1!, A, tolf = options.tolf)
end

"""
    solve_steady_state_core!(context::Context, x0::AbstractVector{T}, J!::Function, A0::AbstractMatrix{T}; tolf = 1e-8) where T <: Real

Call the nonlinear solver to solve for the steady state
"""
function solve_steady_state_core!(
    context::Context,
    x0::AbstractVector{T},
    J!::Function,
    A0::AbstractMatrix{T};
    tolf = 1e-8,
) where {T<:Real}
    ws = StaticWs(context)
    m = context.models[1]
    w = context.work
    params = w.params
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state

    function f!(residuals::AbstractVector{Float64}, x::AbstractVector{Float64})
        DFunctions.static!(ws.temporary_values, residuals, x, exogenous, w.params)
    end

    residuals = zeros(m.endogenous_nbr)
    f!(residuals, x0)
    of = OnceDifferentiable(f!, J!, vec(x0), residuals, A0)
    result = nlsolve(of, x0; method = :robust_trust_region, show_trace = true, ftol = tolf)
    @debug result
    if converged(result)
        results.trends.endogenous_steady_state .= result.zero
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed())
    end
    return nothing
end

"""
    solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}, options)

Solve numerically for the steady state of a Ramsey problem
"""
function solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}, options)
    ws = StaticWs(context)
    m = context.models[1]
    w = context.work
    params = w.params
    endogenous = zeros(m.endogenous_nbr)
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state
    orig_endo_nbr = m.original_endogenous_nbr
    mult_indices =
        orig_endo_nbr .+ findall(a -> a["type"] == 6, context.models[1].auxiliary_variables)
    mult_nbr = length(mult_indices)
    mult = zeros(mult_nbr)
    orig_endo_aux_nbr = mult_indices[1] - 1
    unknown_variable_indices =
        setdiff!(collect(1:m.original_endogenous_nbr), w.analytical_steadystate_variables)
    unknown_variable_nbr = length(unknown_variable_indices)
    M = zeros(orig_endo_nbr, mult_nbr)
    U1 = zeros(orig_endo_nbr)
    residuals = zeros(m.endogenous_nbr)
    x00 = zeros(unknown_variable_nbr)
    x00 .= view(x0, unknown_variable_indices)

    function f!(x::AbstractVector{Float64})
        view(endogenous, unknown_variable_indices) .= x
        # Lagrange multipliers are kept to zero
        context.modfileinfo.has_auxiliary_variables &&
            DFunctions.static_auxiliary_variables!(endogenous, exogenous, params)
        context.modfileinfo.has_steadystate_file &&
            DFunctions.steady_state!(endogenous, exogenous, params)
        residuals = get_static_residuals!(ws, w.params, endogenous, exogenous)
        U1 .= view(residuals, 1:orig_endo_nbr)
        A = get_static_jacobian!(ws, w.params, endogenous, exogenous, m)
        M .= view(A, 1:orig_endo_nbr, mult_indices)
        mult = view(endogenous, mult_indices)
        mult .= -M \ U1
        view(residuals, 1:orig_endo_nbr) .= U1 .+ M * mult
        res1 = sum(x-> x*x, residuals)
        res2 = norm(residuals)
        return res2
    end

    if unknown_variable_nbr == 0
        res = f!(Float64[])
        if res > options.tolf
            @debug "Steady state computation failed"
            throw(DynareSteadyStateComputationFailed())
        else
            results.trends.endogenous_steady_state .= endogenous
        end
    else
        result = optimize(f!, x00, LBFGS(), Optim.Options(f_tol=1e-6))
        @debug result
        @show result
        if Optim.converged(result) && abs(Optim.minimum(result)) < options.tolf
            view(endogenous, unknown_variable_indices) .= Optim.minimizer(result)
            context.modfileinfo.has_auxiliary_variables &&
                context.dynarefunctions.set_auxiliary_variables!(endogenous, exogenous, params)
            context.modfileinfo.has_steadystate_file &&
                context.dynarefunctions.steady_state!(endogenous, exogenous, params)
            # get Lagrance multipliers
            f!(Optim.minimizer(result))
            results.trends.endogenous_steady_state .= endogenous
        else
            @debug "Steady state computation failed with\n $result"
            throw(DynareSteadyStateComputationFailed())
        end
    end
end
