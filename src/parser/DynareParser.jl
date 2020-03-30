using JSON

@enum SymbolType Endogenous Exogenous ExogenousDeterministic Parameter DynareFunction

struct Symbol
    longname::String
    texname::String
    type::SymbolType
    orderintype::Integer
end

struct Options
end

struct Results
end

struct Context
    symboltable::Dict{String, Symbol}
    models::Vector{Model}
    options::Options
    results::Results
    function Context()
        symboltable = Dict()
        models = Vector{Model}(undef,0)
        options = Options()
        results = Results()
        new(symboltable, models, options, results)
    end
end

function parser(modfilename, context::Context)
    modelstring = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)

    endo_nbr = set_symbol_table!(context.symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(context.symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(context.symboltable, modeljson["exogenous_deterministic"], ExogenousDeterministic)
    param_nbr = set_symbol_table!(context.symboltable, modeljson["parameters"], Parameter)
    params = Vector{Float64}(undef, param_nbr)
    Sigma_e = zeros(exo_nbr, exo_nbr)
    for field in modeljson["statements"]
        if field["statementName"] == "param_init"
            initialize_parameter!(params, field, context.symboltable)
        elseif field["statementName"] == "native"
            native_statement(field)
        elseif field["statementName"] == "initval"
            initval(field)
        elseif field["statementName"] == "shocks"
            shocks!(Sigma_e, field, context.symboltable)
        elseif field["statementName"] == "stoch_simul"
            stoch_simul(field)
        elseif field["statementName"] == "verbatim"
            verbatim(field)
        elseif field["statementName"] == "check"
            check(field)
        elseif field["statementName"] == "stoch_simul"
            check(field)
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

function stoch_simul(field)
end

function verbatim(field)
    println("VERBATIM: $field")
end

function check(field)
end
