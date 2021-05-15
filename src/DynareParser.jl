using CSV
using DataFrames
#using .DynareContainers
using FastLapackInterface
using FastLapackInterface.SchurAlgo
using JSON
using KalmanFilterTools
using LinearRationalExpectations
using Periods
using StatsFuns
using TimeDataFrames

context = Context(Dict{String, DynareSymbol}(),
                  Vector{Model}(undef, 0),
                  Dict(),
                  Results(Vector{ModelResults}(undef, 0)),
                  Work(Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Matrix{Float64}(undef, 0, 0),
                       Vector{Int64}(undef, 0),
                       Matrix{Float64}(undef, 0, 0),
                       Matrix{Float64}(undef, 0, 0),
                       false,
                       Matrix{Float64}(undef, 0, 0)
                       )
                  )
                         
function parser(modfilename::String, commandlineoptions::CommandLineOptions)
    modelstring::String = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)
    
    symboltable = SymbolTable()
    endo_nbr = set_symbol_table!(symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(symboltable,
                                    modeljson["exogenous_deterministic"],
                                    ExogenousDeterministic)
    param_nbr = set_symbol_table!(symboltable, modeljson["parameters"], Parameter)
    model_info = get_model_info(modeljson["model_info"])

    model = Model(modfilename,
                  endo_nbr,
                  model_info.lead_lag_incidence,
                  exo_nbr,
                  0,
                  exo_det_nbr,
                  param_nbr,
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
                  model_info.NNZDerivatives,
                  commandlineoptions.compilemodule
                  )
    varobs = Vector{String}()
    if "varobs" in keys(modeljson)
        varobs = vcat(varobs, modeljson["varobs"])
    end
    if "varexobs" in keys(modeljson)
        varobs = vcat(varobs,modeljson["varexobs"])
    end
    if "varexdetobs" in keys(modeljson)
        varobs = vcat(varobs,modeljson["varexdetobs"])
    end

    order = 1
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
                                Dict{String, Any}())
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2*model.n_both
    ncol1 = ncol + model.exogenous_nbr
    nrow = model.maximum_exo_lag + model.maximum_exo_lead + 1
    work = Work(Vector{Float64}(undef, model.parameter_nbr),
                Vector{Float64}(undef, model.endogenous_nbr),
                Vector{Float64}(undef, sum(model.dynamic!.tmp_nbr[1:2])),
                Vector{Float64}(undef, ncol),
                Matrix{Float64}(undef, nrow, model.exogenous_nbr),
                varobs, 
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                false,
                Matrix{Float64}(undef, 0, 0))
    results = Results([modelresults])
    global context = Context(symboltable, [model], Dict(), results, work)
    if "statements" in keys(modeljson)
        parse_statements!(context, modeljson["statements"])
    end
    return context
end

function parse_statements!(context::Context, statements::Vector{Any})
    for field in statements
        if field["statementName"] == "calib_smoother"
            calib_smoother!(context, field)
        elseif field["statementName"] == "check"
            check!(context, field)
        elseif field["statementName"] == "deterministic_trends"
            deterministic_trends!(context, field)
        elseif field["statementName"] == "histval"
            histval!(context, field)
        elseif field["statementName"] == "initval"
            initval!(context, field)
        elseif field["statementName"] == "native"
            dynare_parse_eval(field["string"], context)
        elseif field["statementName"] == "param_init"
            param_init!(context, field)
        elseif field["statementName"] == "perfect_foresight_setup"
            perfect_foresight_setup!(context, field)
        elseif field["statementName"] == "perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
        elseif field["statementName"] == "planner_objective"
            planner_objective!(context, field)
        elseif field["statementName"] == "ramsey_model"
            Nothing
        elseif field["statementName"] == "shocks"
            shocks!(context, field)
        elseif field["statementName"] == "stoch_simul"
            stoch_simul!(context, field)
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
    return eval(e)
end

function dynare_eval(expr::Expr, context::Context)
    for (i, a) in enumerate(expr.args)
        expr.args[i] = dynare_eval(a, context)
    end
    return expr
end

function dynare_eval(s::Symbol, context::Context)
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
                           symboltype::SymbolType)
    count = 0
    for entry in modelfile
        count += 1
        symbol = DynareSymbol(entry["longName"],
                              entry["texName"],
                              symboltype,
                              count)
        table[entry["name"]] = symbol

    end
    return count
end

get_model_info(field::Dict{String, Any}) =
    ModelInfo(hcat(field["lead_lag_incidence"]...),
              field["nstatic"],
              field["nfwrd"],
              field["npred"],
              field["nboth"],
              field["nsfwrd"],
              field["nspred"],
              field["ndynamic"],
              field["maximum_endo_lag"],
              field["maximum_endo_lead"],
              field["maximum_exo_lag"],
              field["maximum_exo_lead"],
              field["maximum_exo_det_lag"],
              field["maximum_exo_det_lead"],
              field["maximum_lag"],
              field["maximum_lead"],
              field["orig_maximum_endo_lag"],
              field["orig_maximum_endo_lead"],
              field["orig_maximum_exo_lag"],
              field["orig_maximum_exo_lead"],
              field["orig_maximum_exo_det_lag"],
              field["orig_maximum_exo_det_lead"],
              max(field["orig_maximum_lag"],
                  field["orig_maximum_lag_with_diffs_expanded"]),
              field["orig_maximum_lead"],
              field["NNZDerivatives"]
              )
                                  
function verbatim(field::Dict{String, Any})
#    println("VERBATIM: $field")
end

function get_smoothed_values(variable_name::String;
                             context::Context=context)
    k = context.symboltable[variable_name].orderintype
    return context.results.model_results[1].smoother["alphah"][k,:]
end
    
