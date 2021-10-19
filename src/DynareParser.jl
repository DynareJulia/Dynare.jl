using CSV
using DataFrames
#using .DynareContainers
using ExtendedDates
using FastLapackInterface
using FastLapackInterface.SchurAlgo
using JSON
using KalmanFilterTools
using LinearRationalExpectations
using Periods
using StatsFuns
using TimeDataFrames

@noinline function parseJSON(modfilename::String)
    modelstring::String = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)
    return modeljson
end

function get_symbol_table(modeljson::Dict{String, Any})
    symboltable = Dict{String, DynareSymbol}()

    endo_nbr = set_symbol_table!(symboltable,
                                 modeljson["endogenous"],
                                 Endogenous)
    exo_nbr = set_symbol_table!(symboltable,
                                modeljson["exogenous"],
                                Exogenous)
    exo_det_nbr = set_symbol_table!(symboltable,
                                    modeljson["exogenous_deterministic"],
                                    ExogenousDeterministic)
    param_nbr = set_symbol_table!(symboltable,
                                  modeljson["parameters"],
                                  Parameter)
    orig_endo_nbr = modeljson["orig_endo_nbr"]::Int64
    aux_vars = modeljson["aux_vars"]
    return (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars)
end

function get_model(modfilename::String,
                   dynare_model_info::Dict{String, Any},
                   commandlineoptions::CommandLineOptions,
                   endo_nbr::Int64, exo_nbr::Int64,
                   exo_det_nbr::Int64, param_nbr::Int64,
                   orig_endo_nbr::Int64, aux_vars::Vector{Any})
    model_info = get_model_info(dynare_model_info)

    NNZDerivatives = Vector{Int64}(undef, length(model_info.NNZDerivatives))
    for (i, n) in enumerate(model_info.NNZDerivatives)
        NNZDerivatives[i] =  n::Int64
    end
    lead_lag_incidence = Vector{Vector{Int64}}(undef, 0)
    for (i, vv) = enumerate(model_info.lead_lag_incidence)
        v = zeros(Int64, 3)
        for j = 1:3
            v[j] = vv[j]::Int64 
        end
        push!(lead_lag_incidence, v)
    end
    model = Model(modfilename,
                  endo_nbr,
                  lead_lag_incidence,
                  exo_nbr,
                  0,
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
                  commandlineoptions.compilemodule
                  )
    return model
end

function get_varobs(modeljson::Dict{String, Any})
    varobs = Vector{String}()
    if "varobs" in keys(modeljson)
        varobs = vcat(varobs, Vector{String}(modeljson["varobs"]))
    end
    if "varexobs" in keys(modeljson)
        varobs = vcat(varobs, Vector{String}(modeljson["varexobs"]))
    end
    if "varexdetobs" in keys(modeljson)
        varobs = vcat(varobs, Vector{String}(modeljson["varexdetobs"]))
    end
    return varobs
end

function make_containers(endo_nbr::Int64, exo_nbr::Int64, exo_det_nbr::Int64, param_nbr::Int64,
                         model::Model, symboltable::SymbolTable, varobs::Vector{String})
    modelresults = ModelResults(Vector{Float64}(undef, endo_nbr),
                                Trends(endo_nbr, exo_nbr, exo_det_nbr),
                                Matrix{Float64}(undef, endo_nbr, endo_nbr),
                                Vector{Bool}(undef, endo_nbr),
                                Vector{Float64}(undef, exo_nbr),
                                Vector{Float64}(undef, exo_det_nbr),
                                LinearRationalExpectationsResults(endo_nbr,
                                                                  exo_nbr,
                                                                  model.n_states),
                                Vector{Simulation}(undef, 0),
                                Dict{String, Matrix{Float64}}())
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2*model.n_both
    ncol1 = ncol + model.exogenous_nbr
    tmp_nbr = model.dynamic!.tmp_nbr::Vector{Int64}
    work = Work(Vector{Float64}(undef, model.parameter_nbr),
                Vector{Float64}(undef, model.endogenous_nbr),
                Vector{Float64}(undef, sum(tmp_nbr[1:2])),
                Vector{Float64}(undef, ncol),
                # reserve enough space for a single period computation
                Vector{Float64}(undef, 3*model.exogenous_nbr),
                varobs, 
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                [false],
                Matrix{Union{Float64, Missing}}(missing,
                                                model.orig_maximum_lag,
                                                model.endogenous_nbr))
    results = Results([modelresults])

    return Context(symboltable, [model], Modelfile(), results, work)
end

function parser(modfilename::String, commandlineoptions::CommandLineOptions)
    modeljson = parseJSON(modfilename)

    (symboltable, endo_nbr, exo_nbr, exo_det_nbr, param_nbr, orig_endo_nbr, aux_vars) =
        get_symbol_table(modeljson)
    model = get_model(modfilename,
                      modeljson["model_info"],
                      commandlineoptions,
                      endo_nbr, exo_nbr, exo_det_nbr,
                      param_nbr, orig_endo_nbr, aux_vars)
    varobs = get_varobs(modeljson)

    global context = make_containers(endo_nbr, exo_nbr, exo_det_nbr,
                                     param_nbr, model, symboltable,
                                     varobs)
    
    if "statements" in keys(modeljson)
        parse_statements!(context, modeljson["statements"])
    end
    return context
end

function parse_statements!(context::Context, statements::Vector{Any})
    modelfile = context.modelfile
    for field in statements
        if field["statementName"] == "calib_smoother"
            calib_smoother!(context, field)
            modelfile["has_calib_smoother"] = true
        elseif field["statementName"] == "check"
            check!(context, field)
            modelfile["has_check"] = true
        elseif field["statementName"] == "deterministic_trends"
            deterministic_trends!(context, field)
            modelfile["has_trends"] = true
        elseif field["statementName"] == "histval"
            histval!(context, field)
            modelfile["has_histval"] = true
        elseif field["statementName"] == "initval"
            initval!(context, field)
            modfile["has_initval"] = true
        elseif field["statementName"] == "native"
            try
                dynare_parse_eval(field["string"], context)
            catch
                error("""Unrecognized statement $(field["statementName"]) $(field["string"])""")
            end
        elseif field["statementName"] == "param_init"
            param_init!(context, field)
        elseif field["statementName"] == "perfect_foresight_setup"
            perfect_foresight_setup!(context, field)
            modelfile["has_perfect_foresight_setup"] = true
        elseif field["statementName"] == "perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
            modelfile["has_perfect_foresight_solver"] = true
        elseif field["statementName"] == "planner_objective"
            planner_objective!(context, field)
            modelfile["has_planner_objective"] = true
        elseif field["statementName"] == "ramsey_model"
            modelfile["has_ramsey_model"] = true
        elseif field["statementName"] == "shocks"
            shocks!(context, field)
            modelfile["has_shocks"] = true
        elseif field["statementName"] == "stoch_simul"
            stoch_simul!(context, field)
            modelfile["has_stoch_simul"] = true
        elseif field["statementName"] == "verbatim"
            Nothing
        else
            error("""Unrecognized statement $(field["statementName"])""")
        end
    end
end

function dynare_parse_eval(s::String, context::Context)
    e = Meta.parse(s)
    e = dynare_eval(e, context)
    try
        return(eval(e))
    catch
        @show s
        rethrow()
    end
end

function dynare_eval(expr::Expr, context::Context)
    for (i, a) in enumerate(expr.args)
        expr.args[i] = dynare_eval(a, context)
    end
    return expr
end

function dynare_eval(s::Symbol, context::Context)::Union{Symbol, Float64}
    symboltable = context.symboltable
    ks = keys(symboltable)
    params = context.work.params
    ss = string(s)
    if ss in ks
        st = symboltable[ss]
        if st.symboltype == Dynare.Parameter
            s = params[st.orderintype]
        end
    end
    return s
end

function dynare_eval(x::Real, context::Context)
    return x
end

function dynare_eval(s::String, context::Context)
    return s
end

function dynare_eval(q::QuoteNode, context::Context)
    return q
end

function set_symbol_table!(table::Dict{String, DynareSymbol},
                           modelfile::Vector{Any},
                           symboltype::SymbolType)::Int64
    count = 0
    for entry in modelfile
        count += 1
        symbol = DynareSymbol(entry["longName"]::String,
                              entry["texName"]::String,
                              symboltype,
                              count)
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

get_model_info(field::Dict{String, Any}) =
    ModelInfo(field["lead_lag_incidence"]::Vector{Any},
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
              max(field["orig_maximum_lag"]::Int64,
                  field["orig_maximum_lag_with_diffs_expanded"]::Int64),
              field["orig_maximum_lead"]::Int64,
              field["NNZDerivatives"]::Vector{Any}
              )
                                  
function verbatim(field::Dict{String, Any})
#    println("VERBATIM: $field")
end

function get_smoothed_values(variable_name::String;
                             context::Context=context)
    k = context.symboltable[variable_name].orderintype
    return context.results.model_results[1].smoother["alphah"][k,:]
end
    
