using CSV
using DataFrames
#using .DynareContainers
using ExtendedDates
using FastLapackInterface
using JLD2
using JSON
using KalmanFilterTools
using LinearRationalExpectations
using Plots
using StatsFuns
using TimeDataFrames

@noinline function parseJSON(modfilename::String)
    modelstring::String =
        open(f -> read(f, String), modfilename * "/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)
    return modeljson
end

function get_symbol_table(modeljson::Dict{String,Any})
    symboltable = Dict{String,DynareSymbol}()

    endo_nbr = set_symbol_table!(symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(
        symboltable,
        modeljson["exogenous_deterministic"],
        ExogenousDeterministic,
    )
    param_nbr = set_symbol_table!(symboltable, modeljson["parameters"], Parameter)
    orig_endo_nbr = modeljson["orig_endo_nbr"]::Int64
    aux_vars = modeljson["aux_vars"]
    return (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars)
end

function get_model(
    modfilename::String,
    modfileinfo::ModFileInfo,
    dynare_model_info::Dict{String,Any},
    commandlineoptions::CommandLineOptions,
    endo_nbr::Int64,
    exo_nbr::Int64,
    exo_det_nbr::Int64,
    param_nbr::Int64,
    orig_endo_nbr::Int64,
    aux_vars::Vector{Any},
    dynamic_g1_sparse_rowval::Vector{Int64},
    dynamic_g1_sparse_colptr::Vector{Int64},
    static_g1_sparse_rowval::Vector{Int64},
    static_g1_sparse_colptr::Vector{Int64},
)
    model_info = get_model_info(dynare_model_info)
    NNZDerivatives = Vector{Int64}(undef, length(model_info.NNZDerivatives))
    for (i, n) in enumerate(model_info.NNZDerivatives)
        NNZDerivatives[i] = n::Int64
    end
    nlags = model_info.maximum_endo_lag + model_info.maximum_endo_lead + 1
    lead_lag_incidence = Vector{Vector{Int64}}(undef, 0)
    for (i, vv) in enumerate(model_info.lead_lag_incidence)
        v = zeros(Int64, nlags)
        for j = 1:nlags
            v[j] = vv[j]::Int64
        end
        push!(lead_lag_incidence, v)
    end
    model = Model(
        modfilename,
        modfileinfo,
        endo_nbr,
        lead_lag_incidence,
        exo_nbr,
        Int64(0),
        exo_det_nbr,
        param_nbr,
        orig_endo_nbr,
        aux_vars,
        model_info.maximum_endo_lag,
        model_info.maximum_endo_lead,
        model_info.maximum_exo_lag,
        model_info.maximum_exo_lead,
        model_info.maximum_exo_det_lag,
        model_info.maximum_exo_det_lead,
        model_info.maximum_lag,
        model_info.maximum_lead,
        model_info.orig_maximum_endo_lag,
        model_info.orig_maximum_endo_lead,
        model_info.orig_maximum_exo_lag,
        model_info.orig_maximum_exo_lead,
        model_info.orig_maximum_exo_det_lag,
        model_info.orig_maximum_exo_det_lead,
        model_info.orig_maximum_lag,
        model_info.orig_maximum_lead,
        NNZDerivatives,
        commandlineoptions.compilemodule,
        dynamic_g1_sparse_rowval,
        dynamic_g1_sparse_colptr,
        static_g1_sparse_rowval,
        static_g1_sparse_colptr
    )
    return model
end

function get_varobs(modeljson::Dict{String,Any})
    varobs = Vector{String}()
    if haskey(modeljson, "varobs")
        varobs = vcat(varobs, Vector{String}(modeljson["varobs"]))
    end
    if haskey(modeljson, "varexobs")
        varobs = vcat(varobs, Vector{String}(modeljson["varexobs"]))
    end
    if haskey(modeljson, "varexdetobs")
        varobs = vcat(varobs, Vector{String}(modeljson["varexdetobs"]))
    end
    return varobs
end

function check_function_files!(modfileinfo::ModFileInfo, modfilename::String)
    if isfile(modfilename * "Dynamic.jl")
        modfileinfo.has_dynamic_file = true
    end
    if isfile(modfilename * "DynamicSetAuxiliarySeries.jl")
        modfileinfo.has_auxiliary_variables = true
        if !isfile(modfilename * "SetAuxiliaryVariables.jl")
            error(modfilename * "SetAuxiliaryVariables.jl is missing")
        end
    end
    if isfile(modfilename * "SteadyState2.jl")
        modfileinfo.has_steadystate_file = true
    end
end

function make_containers(
    modelfileinfo::ModFileInfo,
    modfilename::String,
    endo_nbr::Int64,
    exo_nbr::Int64,
    exo_det_nbr::Int64,
    param_nbr::Int64,
    model::Model,
    symboltable::SymbolTable,
    varobs::Vector{String},
    commandlineoption,
)
    work = Work(model, varobs)
    modelresults = ModelResults(
        Vector{Float64}(undef, endo_nbr),
        Dict{Symbol,TimeDataFrame}(),
        Trends(endo_nbr, exo_nbr, exo_det_nbr),
        Vector{Bool}(undef, endo_nbr),
        Vector{Float64}(undef, exo_nbr),
        Vector{Float64}(undef, exo_det_nbr),
        LinearRationalExpectationsResults(endo_nbr, exo_nbr, model.n_states),
        Vector{Simulation}(undef, 0),
        Dict{String,Matrix{Float64}}(),
    )
    results = Results([modelresults])
    dynarefunctions = DynareFunctions(
        commandlineoption.compilemodule,
        modelfileinfo,
        modfilename,
        model.orig_maximum_lag,
        model.orig_maximum_lead,
    )
    return Context(symboltable, [model], dynarefunctions, modelfileinfo, results, work)
end

function parser(modfilename::String, commandlineoptions::CommandLineOptions)
    @debug "$(now()): Start $(nameof(var"#self#"))"

    modeljson = parseJSON(modfilename)
    @debug "$(now()): get symbol_table"
    (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars) =
        get_symbol_table(modeljson)
    @debug "$(now()): get Modelfile"
    modfileinfo = ModFileInfo(modfilename)
    check_function_files!(modfileinfo, modfilename)
    @debug "$(now()): get model"
    model = get_model(
        modfilename,
        modfileinfo,
        modeljson["model_info"],
        commandlineoptions,
        endo_nbr,
        exo_nbr,
        exo_det_nbr,
        param_nbr,
        orig_endo_nbr,
        aux_vars,
        Vector{Int64}(modeljson["dynamic_g1_sparse_rowval"]),
        Vector{Int64}(modeljson["dynamic_g1_sparse_colptr"]),
        Vector{Int64}(modeljson["static_g1_sparse_rowval"]),
        Vector{Int64}(modeljson["static_g1_sparse_colptr"]),
    )
    varobs = get_varobs(modeljson)
    @debug "$(now()): make_container"
    global context = make_containers(
        modfileinfo,
        modfilename,
        endo_nbr,
        exo_nbr,
        exo_det_nbr,
        param_nbr,
        model,
        symboltable,
        varobs,
        commandlineoptions,
    )
    get_mcps!(context.models[1].mcps, modeljson["model"])
    if haskey(modeljson, "statements")
        parse_statements!(context, modeljson["statements"])
    end
    @info "$(now()): End $(nameof(var"#self#"))"
    return context
end

function parse_statements!(context::Context, statements::Vector{Any})
    @info "$(now()): Start $(nameof(var"#self#"))"
    modfileinfo = context.modfileinfo
    for field in statements
        statementname = field["statementName"]
        if statementname == "calib_smoother"
            calib_smoother!(context, field)
            modfileinfo.has_calib_smoother = true
        elseif statementname == "check"
            check!(context, field)
            modfileinfo.has_check = true
        elseif statementname == "corr_prior"
            parse_prior!(context, field)
        elseif statementname == "deterministic_trends"
            deterministic_trends!(context, field)
            modfileinfo.has_trends = true
        elseif statementname == "estimated_params"
            #            parse_estimated_parameters!(context, field)
            #            modfileinfo.has_trends = true
        elseif statementname == "estimation"
        elseif statementname == "histval"
            histval!(context, field)
            modfileinfo.has_histval = true
        elseif statementname == "initval"
            @debug "$(now()): start initval"
            initval!(context, field)
            modfileinfo.has_initval = true
            @debug "$(now()): end initval"
        elseif statementname == "initval_file"
            @debug "$(now()): start initval_file"
            initval_file!(context, field)
            modfileinfo.has_initval_file = true
            @debug "$(now()): end initval"
        elseif statementname == "native"
            try
                dynare_parse_eval(field["string"], context)
            catch
                error("""Unrecognized statement $(statementname) $(field["string"])""")
            end
        elseif statementname == "param_init"
            param_init!(context, field)
        elseif statementname == "perfect_foresight_setup"
            @debug "$(now()): start perfect_foresight_setup"
            perfect_foresight_setup!(context, field)
            modfileinfo.has_perfect_foresight_setup = true
            @debug "$(now()): end perfect_foresight_setup"
        elseif statementname == "perfect_foresight_solver"
            @debug "$(now()): start perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
            modfileinfo.has_perfect_foresight_solver = true
            @debug "$(now()): end perfect_foresight_solver"
        elseif statementname == "planner_objective"
            planner_objective!(context, field)
            modfileinfo.has_planner_objective = true
        elseif statementname == "prior"
            parse_prior!(context, field)
        elseif statementname == "ramsey_constraints"
            ramsey_constraints!(context, field)
        elseif statementname == "ramsey_model"
            modfileinfo.has_ramsey_model = true
        elseif statementname == "shocks"
            shocks!(context, field)
            modfileinfo.has_shocks = true
        elseif statementname == "std_prior"
            field["name2"] = field["name1"] = field["name"]
            parse_prior!(context, field)
        elseif statementname == "steady"
            @debug "$(now()): start steady"
            steady!(context, field)
            @debug "$(now()): end steady"
        elseif statementname == "stoch_simul"
            @debug "$(now()): start stoch_simul"
            stoch_simul!(context, field)
            modfileinfo.has_stoch_simul = true
            @debug "$(now()): end stoch_simul"
        elseif statementname == "verbatim"
            Nothing
        else
            error("""Unrecognized statement $(statementname)""")
        end
    end
    last_steps(context)
    @info "$(now()): End $(nameof(var"#self#"))"
end

function dynare_parse_eval(
    s::String,
    context::Context;
    xs = SymbolType[],
    xw = Vector{Float64}[],
)
    e = Meta.parse(s)
    e = dynare_eval(e, context, xs, xw)
    try
        return (eval(e))
    catch
        @show s
        rethrow()
    end
end

function dynare_eval(expr::Expr, context::Context, xs, xw)
    for (i, a) in enumerate(expr.args)
        expr.args[i] = dynare_eval(a, context, xs, xw)
    end
    return expr
end

function dynare_eval(s::Symbol, context::Context, xs, xw)::Union{Symbol,Float64}
    symboltable = context.symboltable
    work = context.work
    params = work.params
    ss = string(s)
    v = Union{Symbol,Float64}
    try
        st = symboltable[ss]
        if st.symboltype == Parameter
            v = params[st.orderintype]
        else
            for (xss, xww) in zip(xs, xw)
                if st.symboltype == xss
                    v = xww[st.orderintype]
                end
            end
        end
    catch
        v = s
    end
    return v
end

function dynare_eval(x::Real, context::Context, xs, xw)
    return x
end

function dynare_eval(s::String, context::Context, xs, xw)
    return s
end

function dynare_eval(q::QuoteNode, context::Context, xs, xw)
    return q
end

function set_symbol_table!(
    table::Dict{String,DynareSymbol},
    modelfile::Vector{Any},
    symboltype::SymbolType,
)::Int64
    count = 0
    for entry in modelfile
        count += 1
        symbol = DynareSymbol(
            entry["longName"]::String,
            entry["texName"]::String,
            symboltype,
            count,
        )
        table[entry["name"]::String] = symbol

    end
    return count
end

function get_lead_lag_incidence(ll::Vector{Any})
    n = length(ll)
    m = length(ll[1])
    llm = zeros(Int64, m, n)
    k = 1
    for v::Vector{Int64} in ll
        copyto!(llm, k, v, 1, m)
        k += m
    end
    return llm
end

get_model_info(field::Dict{String,Any}) = ModelInfo(
    field["lead_lag_incidence"]::Vector{Any},
    field["nstatic"]::Int64,
    field["nfwrd"]::Int64,
    field["npred"]::Int64,
    field["nboth"]::Int64,
    field["nsfwrd"]::Int64,
    field["nspred"]::Int64,
    field["ndynamic"]::Int64,
    field["maximum_endo_lag"]::Int64,
    field["maximum_endo_lead"]::Int64,
    field["maximum_exo_lag"]::Int64,
    field["maximum_exo_lead"]::Int64,
    field["maximum_exo_det_lag"]::Int64,
    field["maximum_exo_det_lead"]::Int64,
    field["maximum_lag"]::Int64,
    field["maximum_lead"]::Int64,
    field["orig_maximum_endo_lag"]::Int64,
    field["orig_maximum_endo_lead"]::Int64,
    field["orig_maximum_exo_lag"]::Int64,
    field["orig_maximum_exo_lead"]::Int64,
    field["orig_maximum_exo_det_lag"]::Int64,
    field["orig_maximum_exo_det_lead"]::Int64,
    max(
        field["orig_maximum_lag"]::Int64,
        field["orig_maximum_lag_with_diffs_expanded"]::Int64,
    ),
    field["orig_maximum_lead"]::Int64,
    field["NNZDerivatives"]::Vector{Any},
)

function verbatim(field::Dict{String,Any})
    #    println("VERBATIM: $field")
end

function get_smoothed_values(variable_name::String; context::Context = context)
    k = context.symboltable[variable_name].orderintype
    return context.results.model_results[1].smoother["alphah"][k, :]
end

function display_graphs(filepath::String)
    graphs = joinpath(filepath, "graphs")
    if Sys.islinux()
        if ENV["DESKTOP_SESSION"] == "gnome"
            run(`/usr/bin/eog $graphs`, wait = false)
        elseif ENV["DESKTOP_SESSION"] == "kde"
            run(`/usr/bin/Gwenview $graphs`, wait = false)
        end
    elseif Sys.iswindows()
        for (i, f) in enumerate(readdir(graphs))
            i > 30 && break
            filename = joinpath(graphs, f)
            run(`start $filename`, wait = false)
        end
    elseif Sys.isapple()
        for (i, f) in enumerate(readdir(graphs))
            i > 30 && break
            filename = joinpath(graphs, f)
            run(`open $filename`, wait = false)
        end
    end
end



struct SavedContext
    symboltable::SymbolTable
    models::Vector{Model}
    modfileinfo::ModFileInfo
    results::Results
    work::Work
end

function save_context(context::Context, filepath::String)
    savedcontext = SavedContext(
        context.symboltable,
        context.models,
        context.modfileinfo,
        context.results,
        context.work,
    )
    filename = split(filepath, "/")[end]
    outputpath = mkpath(joinpath(filepath, "output"))
    save(joinpath(outputpath, "$(filename).jld2"), "context", savedcontext)
end


function last_steps(context::Context)
    filepath = context.modfileinfo.modfilepath
    # display graphs
    if isinteractive() &&
       (get(ENV, "TERM_PROGRAM", "") != "vscode") &&
       ("graphs" in readdir(filepath))
        display_graphs(filepath)
    end
    # save context
    save_context(context, filepath)
end
