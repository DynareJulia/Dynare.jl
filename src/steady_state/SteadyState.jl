using Optim

mutable struct DynareSteadyStateComputationFailed <: Exception end
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
    tolf::Float64 - tolerance criterium for equations error [eps()^(1/3)]
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
            elseif k == "tolf"
                tolf = v::Float64
            elseif k == "tolx"
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
    function `steady!`(context::Context, field::Dict{String, Any})

computes the steady state of the model and set the result in `context`
"""
function steady!(context::Context, field::Dict{String,Any})
    modfileinfo = context.modfileinfo
    options = SteadyOptions(get(field, "options", Dict{String,Any}()))
    if (modfileinfo.has_steadystate_file
        && length(context.dynarefunctions.analytical_steady_state_variables) == context.models[1].original_endogenous_nbr)
        @show "compute steady state"
        compute_steady_state!(context)
    else
        results = context.results.model_results[1]
        # will fail if missing values are encountered
        x0 = Float64.(vec(view(context.work.initval_endogenous, 1, :)))
        copy!(results.trends.exogenous_steady_state,
              Float64.(vec(view(context.work.initval_exogenous, 1, :))))
        if modfileinfo.has_ramsey_model
            @show "solve ramsey"
            solve_ramsey_steady_state!(context, x0, options)
        else
            @show "solve steady"
            solve_steady_state!(context, x0, options)
        end
    end
    if options.display
        steadystate_display(context)
    end
end

"""
    function `steadystate_display`(context::Context)

displays the steady state of the model
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
    function `compute_steady_state!`(context::Context)

computes the steady state of the model using solution provided by the user
"""
function compute_steady_state!(context::Context)
    df = context.dynarefunctions
    results = context.results.model_results[1]
    work = context.work
    # explicit steady state
    evaluate_steady_state!(results, df.steady_state!, work.params)
end

"""
    function evaluate_steady_state!(results::ModelResults,
                                static_module::Module,
                                params::AbstractVector{Float64})

evaluates the steady state function provided by the user
"""
function evaluate_steady_state!(
    results::ModelResults,
    steady_state!::Function,
    params::AbstractVector{Float64},
)
    fill!(results.trends.exogenous_steady_state, 0.0)
    steady_state!(results.trends.endogenous_steady_state,
                  results.trends.exogenous_steady_state,
                  params)
end

function sparse_static_jacobian(ws, params, x, exogenous, m, df)
    J = get_static_jacobian!(ws, params, x, exogenous, m, df)
    return sparse(J)
end


"""
    function solve_steady_state!(context::Context,
                                 x0::Vector{Float64})

solves the static model to obtain the steady state
"""
function solve_steady_state!(context::Context, x0::AbstractVector{Float64}, options)

    ws = StaticWs(context)
    m = context.models[1]
    df = context.dynarefunctions
    w = context.work
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state

    if count(!iszero, get_static_jacobian!(ws, w.params, x0, exogenous, m, df)) <
       0.1 * m.endogenous_nbr * m.endogenous_nbr
        function J1!(A, x::AbstractVector{Float64})
            A .= sparse_static_jacobian(ws, w.params, x, exogenous, m, df)
        end
        A0 = sparse_static_jacobian(ws, w.params, x0, exogenous, m, df)
        solve_steady_state_core!(context, x0, J1!, A0, tolf = options.tolf)
    else
        function J2!(A::AbstractMatrix{Float64}, x::AbstractVector{Float64})
            A .= get_static_jacobian!(ws, w.params, x, exogenous, m, df)
        end
        A0 = Matrix{Float64}(undef, m.endogenous_nbr, m.endogenous_nbr)
        J2!(A0, x0)
        solve_steady_state_core!(context, x0, J2!, A0, tolf = options.tolf)
    end
end

function solve_steady_state_core!(context, x0, J!, A0; tolf = 1e-8)

    ws = StaticWs(context)
    m = context.models[1]
    df = context.dynarefunctions
    w = context.work
    params = w.params
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state

    function f!(residuals::AbstractVector{Float64}, x::AbstractVector{Float64})
        context.modfileinfo.has_auxiliary_variables &&
            context.dynarefunctions.set_auxiliary_variables!(x, exogenous, params)
        residuals .= get_static_residuals!(ws, params, x, exogenous, df)
    end

    residuals = zeros(m.endogenous_nbr)
    f!(residuals, x0)
    of = OnceDifferentiable(f!, J!, vec(x0), residuals, A0)
    result = nlsolve(of, x0; method = :robust_trust_region, show_trace = true, ftol = tolf)
    if converged(result)
        results.trends.endogenous_steady_state .= result.zero
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed)
    end
end

function solve_ramsey_steady_state!(context::Context, x0::AbstractVector{Float64}, options)
    ws = StaticWs(context)
    m = context.models[1]
    df = context.dynarefunctions
    w = context.work
    params = w.params
    endogenous = zeros(m.endogenous_nbr)
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state
    orig_endo_nbr = m.original_endogenous_nbr
    mult_indices = orig_endo_nbr .+ findall(a -> a["type"] == 6, context.models[1].auxiliary_variables)
    mult_nbr = length(mult_indices)
    mult = zeros(mult_nbr)
    orig_endo_aux_nbr = mult_indices[1] - 1
    unknown_variable_indices = setdiff!(collect(1:m.original_endogenous_nbr), df.analytical_steady_state_variables)
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
            context.dynarefunctions.set_auxiliary_variables!(endogenous, exogenous, params)
        context.modfileinfo.has_steadystate_file &&
            context.dynarefunctions.steady_state!(endogenous, exogenous, params)
        residuals = get_static_residuals!(ws, w.params, endogenous, exogenous, df)
        U1 .= view(residuals, 1:orig_endo_nbr)
        A = get_static_jacobian!(ws, w.params, endogenous, exogenous, m, df)
        M .= view(A, 1:orig_endo_nbr, mult_indices)
        mult = -M\U1
        view(residuals, 1:orig_endo_nbr) .= U1 .+ M*mult
        return sum(z -> z*z, residuals)
    end

    result = optimize(f!, x00, LBFGS())
    @show result
    if Optim.converged(result) && abs(result.minimum) < 1e-8
        view(endogenous, unknown_variable_indices) .= result.minimizer
        context.modfileinfo.has_auxiliary_variables &&
            context.dynarefunctions.set_auxiliary_variables!(endogenous, exogenous, params)
        context.modfileinfo.has_steadystate_file &&
            context.dynarefunctions.steady_state!(endogenous, exogenous, params)
        results.trends.endogenous_steady_state .= endogenous
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed)
    end
end
