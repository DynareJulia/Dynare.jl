using CSV
using DataFrames
using ExtendedDates
using FastLapackInterface
using JSON
using Plots
using StatsFuns
using AxisArrayTables

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
    orig_endo_nbr = modeljson["orig_endo_nbr"]::Int
    aux_vars = modeljson["aux_vars"]
    return (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars)
end

function get_model(
    modfilename::String,
    modfileinfo::ModFileInfo,
    dynare_model_info::Dict{String,Any},
    commandlineoptions::CommandLineOptions,
    endo_nbr::Int,
    exo_nbr::Int,
    exo_det_nbr::Int,
    param_nbr::Int,
    orig_endo_nbr::Int,
    aux_vars::Vector{Any},
    dynamic_g1_sparse_rowval::Vector{Int},
    dynamic_g1_sparse_colval::Vector{Int},
    dynamic_g1_sparse_colptr::Vector{Int},
    dynamic_g2_sparse_indices::Vector{Vector{Int}},
    static_g1_sparse_rowval::Vector{Int},
    static_g1_sparse_colptr::Vector{Int},
    dynamic_tmp_nbr::Vector{Int},
    static_tmp_nbr::Vector{Int}
)
    model_info = get_model_info(dynare_model_info)
    NNZDerivatives = Vector{Int}(undef, length(model_info.NNZDerivatives))
    for (i, n) in enumerate(model_info.NNZDerivatives)
        NNZDerivatives[i] = n::Int
    end
    nlags = model_info.maximum_endo_lag + model_info.maximum_endo_lead + 1
    lead_lag_incidence = Vector{Vector{Int}}(undef, 0)
    for (i, vv) in enumerate(model_info.lead_lag_incidence)
        v = zeros(Int, nlags)
        for j = 1:nlags
            v[j] = vv[j]::Int
        end
        push!(lead_lag_incidence, v)
    end
    model = Model(
        modfilename,
        modfileinfo,
        endo_nbr,
        lead_lag_incidence,
        exo_nbr,
        Int(0),
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
        dynamic_g1_sparse_colval,
        dynamic_g1_sparse_colptr,
        dynamic_g2_sparse_indices,
        static_g1_sparse_rowval,
        static_g1_sparse_colptr,
        dynamic_tmp_nbr,
        static_tmp_nbr
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

function isfile_recent(filename, modfilename)
    pathname = joinpath(modfilename, "model/julia", filename)
    return isfile(pathname)
end 

function check_function_files!(modfileinfo::ModFileInfo, modfilename::String)
    dirname = modfilename * "/model/julia/"
    modfileinfo.has_dynamic_file = isfile_recent("SparseDynamicResid!.jl", modfilename)
    modfileinfo.has_static_file = isfile_recent("SparseStaticResid!.jl", modfilename)
    modfileinfo.has_steadystate_file = isfile_recent("SteadyState2.jl", modfilename)
    modfileinfo.has_auxiliary_variables = isfile_recent("DynamicSetAuxiliarySeries.jl", modfilename)
    if modfileinfo.has_auxiliary_variables && !isfile(dirname * "SetAuxiliaryVariables.jl")
        error(dirname * "SetAuxiliaryVariables.jl is missing")
    end
end

function make_containers(
    modelfileinfo::ModFileInfo,
    modfilename::String,
    endo_nbr::Int,
    exo_nbr::Int,
    exo_det_nbr::Int,
    param_nbr::Int,
    model::Model,
    symboltable::SymbolTable,
    varobs::Vector{String},
    commandlineoption,
)
    work = Work(model, varobs)
    modelresults = ModelResults(
        Dict{Symbol, AxisArrayTable}(),
        Trends(endo_nbr, exo_nbr, exo_det_nbr),
        Vector{Bool}(undef, endo_nbr),
        EstimationResults(),
        AxisArrayTable(Matrix{Float64}(undef, 0, 0), [], Symbol[]),
        [AxisArrayTable(Matrix{Float64}(undef, 0, 0), [], Symbol[])],
        AxisArrayTable(Matrix{Float64}(undef, 0, 0), [], Symbol[]),
        LinearRationalExpectationsResults(endo_nbr, exo_nbr, model.n_states),
        Vector{Simulation}(undef, 0),
        AxisArrayTable(Matrix{Float64}(undef, 0, 0), [], Symbol[]),
        Vector{Matrix{Float64}}(undef, 0)
    )
    results = Results([modelresults])
    return Context(symboltable, [model], modelfileinfo, results, work, Dict())
end

function make_context(modeljson, modfilename, commandlineoptions)
    @debug "$(now()): get symbol_table"
    (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars) =
        get_symbol_table(modeljson)
    @debug "$(now()): get Modelfile"
    modfileinfo = ModFileInfo(modfilename)
    check_function_files!(modfileinfo, modfilename)
    modfileinfo.has_steadystate_file = modeljson["steady_state_model"]
    @debug "$(now()): get model"
    dynamic_g2_sparse_indices = Vector{Vector{Int}}(get(modeljson,
                                                        "dynamic_g2_sparse_indices",
                                                        [])
                                                    )
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
        Vector{Int}(modeljson["dynamic_g1_sparse_rowval"]),
        Vector{Int}(modeljson["dynamic_g1_sparse_colval"]),
        Vector{Int}(modeljson["dynamic_g1_sparse_colptr"]),
        dynamic_g2_sparse_indices,
        Vector{Int}(modeljson["static_g1_sparse_rowval"]),
        Vector{Int}(modeljson["static_g1_sparse_colptr"]),
        Vector{Int}(modeljson["dynamic_tmp_nbr"]),
        Vector{Int}(modeljson["static_tmp_nbr"]),
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
    return context
end

function parser(modfilename::String, commandlineoptions::CommandLineOptions)
    @debug "$(now()): Start $(nameof(var"#self#"))"

    modeljson = parseJSON(modfilename)
    context = make_context(modeljson, modfilename, commandlineoptions)
    context.work.analytical_steadystate_variables = DFunctions.load_model_functions(modfilename)
    if haskey(modeljson, "statements")
        if commandlineoptions.stoponerror
            parse_statements!(context, modeljson["statements"])
        else    
            try
                parse_statements!(context, modeljson["statements"])
            catch e
                println(e)
            end
        end 
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
        elseif statementname == "deterministic_trends"
            deterministic_trends!(context, field)
            modfileinfo.has_trends = true
        elseif statementname == "endval"
            @debug "$(now()): start endval"
            endval!(context, field)
            modfileinfo.has_endval = true
            modfileinfo.endval_is_reset = true
            @debug "$(now()): end endval"
        elseif statementname == "estimated_params"
            parse_estimated_parameters!(context, field)
        elseif statementname == "estimated_params_init"
            parse_estimated_parameter_init!(context, field)
        elseif statementname == "estimation"
            estimation!(context, field)
        elseif statementname == "histval"
            histval!(context, field)
            modfileinfo.has_histval = true
        elseif statementname == "homotopy"
            homotopy_setup!(context, field)
        elseif statementname == "initval"
            @debug "$(now()): start initval"
            initval!(context, field)
            modfileinfo.has_initval = true
            modfileinfo.initval_is_reset = true
            @debug "$(now()): end initval"
        elseif statementname == "initval_file"
            @debug "$(now()): start initval_file"
            initval_file!(context, field)
            modfileinfo.has_initval_file = true
            @debug "$(now()): end initval"
        elseif statementname == "load_params_and_steady_state"
            load_params_and_steady_state(context, field)
        elseif statementname == "native"
            try
                dynare_parse_eval(field["string"], context)
            catch
                error("""Unrecognized statement $(statementname) $(field["string"])""")
            end
        elseif statementname == "param_init"
            # implicit parameter initialization in the *.mod file
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
        elseif statementname == "ramsey_constraints"
            ramsey_constraints!(context, field)
        elseif statementname == "ramsey_model"
            modfileinfo.has_ramsey_model = true
        elseif statementname == "save_params_and_steady_state"
            save_params_and_steady_state(context, field)
        elseif statementname == "shocks"
            shocks!(context, field)
            modfileinfo.has_shocks = true
        elseif statementname == "steady"
            @debug "$(now()): start steady"
            steady_!(context, field)
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
    results = context.results.model_results[1]
    trends = results.trends
    ss = string(s)
    let v::Union{Symbol,Float64}
        try
            st = symboltable[ss]
            k = st.orderintype
            if st.symboltype == Parameter
                v = params[k]
            elseif st.symboltype == Exogenous
                v = results.exogenous_steady_state[k]
            else
                for (xss, xww) in zip(xs, xw)
                    if st.symboltype == xss
                        v = xww[k]
                    end
                end
            end
        catch
            v = s
        end
        return v
    end 
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
)::Int
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
    llm = zeros(Int, m, n)
    k = 1
    for v::Vector{Int} in ll
        copyto!(llm, k, v, 1, m)
        k += m
    end
    return llm
end

get_model_info(field::Dict{String,Any}) = ModelInfo(
    field["lead_lag_incidence"]::Vector{Any},
    field["nstatic"]::Int,
    field["nfwrd"]::Int,
    field["npred"]::Int,
    field["nboth"]::Int,
    field["nsfwrd"]::Int,
    field["nspred"]::Int,
    field["ndynamic"]::Int,
    field["maximum_endo_lag"]::Int,
    field["maximum_endo_lead"]::Int,
    field["maximum_exo_lag"]::Int,
    field["maximum_exo_lead"]::Int,
    field["maximum_exo_det_lag"]::Int,
    field["maximum_exo_det_lead"]::Int,
    field["maximum_lag"]::Int,
    field["maximum_lead"]::Int,
    field["orig_maximum_endo_lag"]::Int,
    field["orig_maximum_endo_lead"]::Int,
    field["orig_maximum_exo_lag"]::Int,
    field["orig_maximum_exo_lead"]::Int,
    field["orig_maximum_exo_det_lag"]::Int,
    field["orig_maximum_exo_det_lead"]::Int,
    max(
        field["orig_maximum_lag"]::Int,
        field["orig_maximum_lag_with_diffs_expanded"]::Int,
    ),
    field["orig_maximum_lead"]::Int,
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
            run(`cmd /k start $filename`, wait = false)
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
    workspaces::Dict
end

function save_context(context::Context, filepath::String)
    filename = split(filepath, "/")[end]
    outputpath = mkpath(joinpath(filepath, "output"))
    serialize(joinpath(outputpath, "$(filename).jls"), context)
    return nothing
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
    return context
end

function get_mcps!(mcps::Vector{Tuple{Int,Int,String,String}},
                   model::Vector{Any})
    for (i, eq) in enumerate(model)
        tags = get(eq, "tags", "")
        if !isempty(tags)
            tag = get(eq["tags"], "mcp", "")
            if !isempty(tag)
                p1, p2, p3 = split(tag, limit = 3)
                iv = context.symboltable[p1].orderintype
                push!(mcps, (i, iv, p2, p3))
            end
        end
    end
end

function ramsey_constraints!(context, field)
    for c in field["ramsey_model_constraints"]
        p1, p2, p3 = split(c["constraint"], limit = 3)
        eq = context.symboltable[p1].orderintype
        push!(context.models[1].mcps, (eq, eq, p2, p3))
    end
end

# redefinition normcdf, normpdf
normcdf(x; μ = 0, σ = 1) = cdf(Normal(μ, σ), x)
normpdf(x; μ = 0, σ = 1) = pdf(Normal(μ, σ), x)
