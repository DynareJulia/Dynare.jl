using FastLapackInterface
using FastLapackInterface.LinSolveAlgo
using JSON
using LinearRationalExpectations

@enum SymbolType Endogenous Exogenous ExogenousDeterministic Parameter DynareFunction

struct Symbol
    longname::String
    texname::String
    type::SymbolType
    orderintype::Integer
end

struct Options
end

struct ModelResults
    endogenous_steady_state::Vector{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_deterministic_steady_state::Vector{Float64}
end

struct Results
    model_results::Vector{ModelResults}
end

mutable struct Work
    params::Vector{Float64}
    residuals::Vector{Float64}
    temporary_values::Vector{Float64}
    dynamic_variables::Vector{Float64}
    jacobian::Matrix{Float64}
    qr_jacobian::Matrix{Float64}
end

function work_allocate(w::Work, m::Model)
    if length(w.params) == 0
        resize!(w.params, m.parameter_nbr)
    end
    if length(w.residuals) == 0
        resize!(w.residuals, m.endogenous_nbr)
    end
    if length(w.temporary_values) == 0
        resize!(w.temporary_values, sum(m.dynamic!.tmp_nbr))
    end
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both
    if length(w.dynamic_variables) == 0
        resize!(w.dynamic_variables, ncol)
    end
    ncol1 = ncol + m.exogenous_nbr
    if length(w.jacobian) == 0
        w.jacobian = Matrix{Float64}(undef, m.endogenous_nbr, ncol1)
    end
    if length(w.qr_jacobian) == 0
        w.qr_jacobian = Matrix{Float64}(undef, m.endogenous_nbr, ncol1)
    end
end          

mutable struct Context
    symboltable::Dict{String, Symbol}
    models::Vector{Model}
    options::Dict{String, Any}
    results::Results
    work::Work
    function Context()
        symboltable = Dict()
        models = Vector{Model}(undef,1)
        options = Dict()
        results = Results([ModelResults([],[],[])])
        work = Work([], [], [], [], Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0))
        new(symboltable, models, options, results, work)
    end
end

struct ModelInfo
    lead_lag_incidence::Array{Int64}
    nstatic
    nfwrd
    npred
    nboth
    nsfwrd
    nspred
    ndynamic
end

function parser(modfilename, context::Context)
    modelstring = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)

    endo_nbr = set_symbol_table!(context.symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(context.symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(context.symboltable, modeljson["exogenous_deterministic"], ExogenousDeterministic)
    param_nbr = set_symbol_table!(context.symboltable, modeljson["parameters"], Parameter)
    Sigma_e = zeros(exo_nbr, exo_nbr)
    model_info = get_model_info(modeljson["model_info"])
    context.models[1] = Model(modfilename,
                              endo_nbr,
                              model_info.lead_lag_incidence,
                              exo_nbr,
                              0,
                              exo_det_nbr,
                              param_nbr)
    context.results = Results([ModelResults(Vector{Float64}(undef, endo_nbr),
                                    Vector{Float64}(undef, exo_nbr),
                                            Vector{Float64}(undef, exo_det_nbr))])
    work_allocate(context.work, context.models[1])
    for field in modeljson["statements"]
        if field["statementName"] == "param_init"
            initialize_parameter!(context.work.params, field, context.symboltable)
        elseif field["statementName"] == "native"
            native_statement(field)
        elseif field["statementName"] == "initval"
            initval(field)
        elseif field["statementName"] == "shocks"
            shocks!(Sigma_e, field, context.symboltable)
        elseif field["statementName"] == "verbatim"
            verbatim(field)
        elseif field["statementName"] == "check"
            check(field)
        elseif field["statementName"] == "stoch_simul"
            println(field)
            stoch_simul!(context, field)
        elseif field["statementName"] == "perfect_foresight_setup"
            perfect_foresight_setup!(options, field)
        elseif field["statementName"] == "perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
        else
            error("Unrecognized statement $(field["statementName"])")
        end
    end
end

function set_symbol_table!(table::Dict{String, Dynare.Symbol},
                     modelfile,
                     type::SymbolType)
    count = 0
    for entry in modelfile
        count += 1
        symbol = Dynare.Symbol(entry["longName"],
                        entry["texName"],
                        type,
                        count)
        table[entry["name"]] = symbol

    end
    return count
end

get_model_info(field) = ModelInfo(hcat(field["lead_lag_incidence"]...),
                                  field["nstatic"],
                                  field["nfwrd"],
                                  field["npred"],
                                  field["nboth"],
                                  field["nsfwrd"],
                                  field["nspred"],
                                  field["ndynamic"])
                   
function initialize_parameter!(params, field, symboltable)
    s = symboltable[field["name"]]
    k = s.orderintype
    params[k] = eval(Meta.parse(field["value"]))
end

function native_statement(field)
    println("NATIVE: $field")
end

function initval(field)
end

function shocks!(Sigma, field, symboltable)
    set_variance!(Sigma, field["variance"], symboltable)
    set_stderr!(Sigma, field["stderr"], symboltable)
    set_covariance!(Sigma, field["covariance"], symboltable)
    set_correlation!(Sigma, field["correlation"], symboltable)
end

function set_variance!(Sigma, variance, symboltable)
    for v in variance
        k =  symboltable[v["name"]].orderintype
        Sigma[k, k] = eval(Meta.parse(v["variance"]))
    end
end

function set_stderr!(Sigma, stderr, symboltable)
    for s in stderr
        k =  symboltable[s["name"]].orderintype
        Sigma[k, k] = eval(Meta.parse(s["stderr"]))
    end
end

function set_covariance!(Sigma, covariance, symboltable)
    for c in covariance
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        Sigma[k1, k2] = eval(Meta.parse(c["covariance"]))
        Sigma[k2, k1] = Sigma[k1, k2]
    end
end

function set_correlation!(Sigma, correlation, symboltable)
    for c in correlation
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        corr = eval(Meta.parse(c["correlation"]))
        Sigma[k2, k1] = sqrt(Sigma[k1, k1]*Sigma[k2, k2])*corr
    end
end

function stoch_simul!(context, field)
    context.options["stoch_simul"] = Dict()
    copy!(context.options["stoch_simul"], field["options"])
    compute_stoch_simul(context)
end

function perfect_foresight_setup(context, field)
    context.options["perfect_foresight_setup"] = Dict()
    copy!(context.options["perfect_foresight_setup"], field["options"])
    compute_perfect_foresight_setup(context)
end

function perfect_foresight_solver(context, field)
    context.options["perfect_foresight_solver"] = Dict()
    copy!(context.options["perfect_foresight_solver"], field["options"])
    compute_perfect_foresight_solver(context)
end


function verbatim(field)
    println("VERBATIM: $field")
end

function check(field)
end

function compute_stoch_simul(context)
    m = context.models[1]
    results = context.results
    if context.options["stoch_simul"]["dr_cycle_reduction"]
        algo = "CR"
    else
        algo = "GS"
    end
    ws = LinearRationalExpectationsWs(algo,
                                      m.endogenous_nbr,
                                      m.exogenous_nbr,
                                      m.exogenous_deterministic_nbr,
                                      m.i_fwrd_b,
                                      m.i_current,
                                      m.i_bkwrd_b,
                                      m.i_both,
                                      m.i_static)
    LinearRationalExpectations.remove_static!(jacobian, ws)
    if algo == "GS"
        LinearRationalExpectations.get_de!(ws, jacobian)
    else
        LinearRationalExpectations.get_abc!(ws, jacobian)
    end
    LinearRationalExpectations.first_order_solver!(results, algo, jacobian, options, ws)
    println(results)
end

function compute_prefect_foresight_setup(context); end
function compute_perfect_foresight_solver(context); end

function get_dynamic_variables!(y::Vector{Float64}, steadystate::Vector{Float64}, lli::Matrix{Int64})
    for i = 1:size(lli,2)
        value = steadystate[i]
        for j = 1:size(lli,1)
            k = lli[j, i]
            if k > 0
                y[k] = value
            end
        end
    end
end

function get_jacobian_at_steadystate!(work::Work, steadystate, exogenous, m::Model, period::Int64)
    lli = m.lead_lag_incidence
    get_dynamic_variables!(work.dynamic_variables, steadystate, lli)
    m.dynamic!.dynamic!(work.temporary_values,
                        work.residuals,
                        work.jacobian,
                        work.dynamic_variables,
                        exogenous,
                        work.params,
                        steadystate,
                        period)  
end
