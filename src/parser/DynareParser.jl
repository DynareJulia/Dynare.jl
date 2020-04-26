using JSON
using LinearRationalExpectations

struct ModelResults
    endogenous_steady_state::Vector{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_deterministic_steady_state::Vector{Float64}
    linearrationalexpectations::LinearRationalExpectationsResults
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

mutable struct Context
    symboltable::Dict{String, DynareSymbol}
    models::Vector{Model}
    options::Dict
    results::Results
    work::Work
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

function parser(modfilename)
    modelstring = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)

    symboltable = SymbolTable()
    endo_nbr = set_symbol_table!(symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(symboltable,
                                    modeljson["exogenous_deterministic"],
                                    ExogenousDeterministic)
    param_nbr = set_symbol_table!(symboltable, modeljson["parameters"], Parameter)
    Sigma_e = zeros(exo_nbr, exo_nbr)
    model_info = get_model_info(modeljson["model_info"])
    model = Model(modfilename,
                  endo_nbr,
                  model_info.lead_lag_incidence,
                  exo_nbr,
                  0,
                  exo_det_nbr,
                  param_nbr)
    order = 1
    modelresults = ModelResults(Vector{Float64}(undef, endo_nbr),
                                Vector{Float64}(undef, exo_nbr),
                                Vector{Float64}(undef, exo_det_nbr),
                                LinearRationalExpectationsResults(order,
                                                                  endo_nbr,
                                                                  exo_nbr,
                                                                  model.n_states))
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2*model.n_both
    ncol1 = ncol + model.exogenous_nbr
    work = Work(Vector{Float64}(undef, model.parameter_nbr),
                Vector{Float64}(undef, model.endogenous_nbr),
                Vector{Float64}(undef, sum(model.dynamic!.tmp_nbr[1:2])),
                Vector{Float64}(undef, ncol),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1))
    results = Results([modelresults])
    context = Context(symboltable, [model], Dict(), results, work)
    for field in modeljson["statements"]
        if field["statementName"] == "param_init"
            initialize_parameter!(work.params, field, symboltable)
        elseif field["statementName"] == "native"
            native_statement(field)
        elseif field["statementName"] == "initval"
            initval(field)
        elseif field["statementName"] == "shocks"
            shocks!(Sigma_e, field, symboltable)
        elseif field["statementName"] == "verbatim"
            verbatim(field)
        elseif field["statementName"] == "check"
            check(field)
        elseif field["statementName"] == "stoch_simul"
            stoch_simul!(context, field)
        elseif field["statementName"] == "perfect_foresight_setup"
            perfect_foresight_setup!(options, field)
        elseif field["statementName"] == "perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
        else
            error("Unrecognized statement $(field["statementName"])")
        end
    end
    return context
end

function set_symbol_table!(table::Dict{String, DynareSymbol},
                     modelfile,
                     type::SymbolType)
    count = 0
    for entry in modelfile
        count += 1
        symbol = DynareSymbol(entry["longName"],
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
#    println("NATIVE: $field")
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

function display_stoch_simul(x, title, context)
    endogenous_names = get_endogenous_longname(context.symboltable)
    emptyrow = ["" for _= 1:size(x,1)]
    column_header = []
    map(x -> push!(column_header, string(x, "\U0209C")), endogenous_names)
    row_header = [""]
    map(x -> push!(row_header, "ϕ($x)"), endogenous_names[context.models[1].i_bkwrd_b])
    map(x -> push!(row_header, "$x\U0209C"), get_exogenous_longname(context.symboltable))
    data = hcat(row_header,
                vcat(reshape(column_header, 1, length(column_header)),
                     x))
    # Note: ϕ(x) = x_{t-1} - \bar x
    note = string("Note: ϕ(x) = x\U0209C\U0208B\U02081 - ", "\U00305", "x")
    println("\n")
    dynare_table(data, title, column_header, row_header, note)
end

function stoch_simul!(context, field)
    context.options["stoch_simul"] = Dict()
    copy!(context.options["stoch_simul"], field["options"])
    compute_stoch_simul(context)
    x = context.results.model_results[1].linearrationalexpectations.g1
    vx = view(x, :, 1:size(x, 2) - 1)
    display_stoch_simul(vx', "Coefficients of approximate solution function", context)
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
#    println("VERBATIM: $field")
end

function check(field)
end

function compute_stoch_simul(context)
    m = context.models[1]
    results = context.results.model_results[1]
    options = context.options["stoch_simul"]
    options["cyclic_reduction"] = Dict()
    options["generalized_schur"] = Dict()
    work = context.work
    Base.invokelatest(steady_state!, context)
    fill!(results.exogenous_steady_state, 0.0)
    get_jacobian_at_steadystate!(work,
                                 results.endogenous_steady_state,
                                 results.exogenous_steady_state, 
                                 m,
                                 2)
    if options["dr_cycle_reduction"]
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
    LinearRationalExpectations.remove_static!(work.jacobian, ws)
    if algo == "GS"
        LinearRationalExpectations.get_de!(ws, work.jacobian)
        options["generalized_schur"]["criterium"] = 1 + 1e-6
    else
        LinearRationalExpectations.get_abc!(ws, work.jacobian)
        options["cyclic_reduction"]["tol"] = 1e-8
    end
    LinearRationalExpectations.first_order_solver!(results.linearrationalexpecations,
                                                   algo,
                                                   work.jacobian,
                                                   options,
                                                   ws)
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
    Base.invokelatest(m.dynamic!.dynamic!,
                      work.temporary_values,
                      work.residuals,
                      work.jacobian,
                      work.dynamic_variables,
                      repeat(exogenous', 3, 1),
                      work.params,
                      steadystate,
                      period)  
end
