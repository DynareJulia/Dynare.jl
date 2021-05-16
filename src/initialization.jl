function histval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    histval = zeros(m.orig_maximum_lag, m.endogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]].orderintype
        l = m.orig_maximum_lag - v["lag"]
        histval[l, k] = dynare_parse_eval(v["value"], context)
    end
    context.work.histval .= histval
end

function param_init!(context::Context, field::Dict{String, Any})
    params = context.work.params;
    symboltable = context.symboltable
    s = symboltable[field["name"]]
    k = s.orderintype
    params[k] = dynare_parse_eval(field["value"], context)
end

function initval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    endo_steady_state = zeros(m.endogenous_nbr)
    exo_steady_state = zeros(m.exogenous_nbr)
    exo_det_steady_state = zeros(m.exogenous_deterministic_nbr)
    for v in field["vals"]
        s = symboltable[v["name"]]
        typ = s.symboltype
        k = s.orderintype
        if typ == :Endogenous
            endo_steady_state[k] = dynare_parse_eval(v["value"], context)
        elseif typ == :Exogenous
            exo_steady_state[k] = dynare_parse_eval(v["value"], context)
        elseif typ == :ExogenousDeterministic
            exo_det_steady_state[k] = dynare_parse_eval(v["value"], context)
        else
            throw(error("$(v["name"]) can't be set in INITVAL"))
        end
    end
    params = context.work.params
    if !isempty(m.set_auxiliary_variables)
        Base.invokelatest(
            m.set_auxiliary_variables!.set_auxiliary_variables!(
                endogenous_steady_state,
                exogenous_steady_state,
                params)
        )
    end
    trends = context.results.results_model[1].trends 
    trends.endogenous_steady_state = endo_steady_state
    trends.exogenous_steady_state = exo_steady_state
    trends.exogenous_det_steady_state = exo_det_steady_state
end

function shocks!(context::Context, field::Dict{String, Any})
    Sigma = context.models[1].Sigma_e
    symboltable = context.symboltable
    set_variance!(Sigma, field["variance"], symboltable)
    set_stderr!(Sigma, field["stderr"], symboltable)
    set_covariance!(Sigma, field["covariance"], symboltable)
    set_correlation!(Sigma, field["correlation"], symboltable)
end

function set_variance!(Sigma::Matrix{Float64}, variance::Vector{Any}, symboltable::SymbolTable)
    for v in variance
        k =  symboltable[v["name"]].orderintype
        Sigma[k, k] = dynare_parse_eval(v["variance"], context)
    end
end

function set_stderr!(Sigma::Matrix{Float64}, stderr::Vector{Any}, symboltable::SymbolTable)
    for s in stderr
        k =  symboltable[s["name"]].orderintype
        x = dynare_parse_eval(s["stderr"], context)
        Sigma[k, k] = x*x
    end
end

function set_covariance!(Sigma::Matrix{Float64}, covariance::Vector{Any}, symboltable::SymbolTable)
    for c in covariance
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        Sigma[k1, k2] = dynare_parse_eval(c["covariance"], context)
        Sigma[k2, k1] = Sigma[k1, k2]
    end
end

function set_correlation!(Sigma::Matrix{Float64}, correlation::Vector{Any}, symboltable::SymbolTable)
    for c in correlation
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        corr = dynare_parse_eval(c["correlation"], context)
        Sigma[k2, k1] = sqrt(Sigma[k1, k1]*Sigma[k2, k2])*corr
    end
end

function load_params!(context::Context, filename::String)
    param_nbr = context.models[1].parameter_nbr;
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = [];
        param_names = get
        for line in readlines(io)
            elem = split(line);
            if elem[1] in keys(symboltable)
                s = symboltable[elem[1]]
                if s.symboltype == Parameter
                    parameters[s.orderintype] = parse(Float64, elem[2])
                end
            end
        end
    end
end

function load_steadystate!(context::Context, filename::String)
    endogenous = context.results.model_results[1].trends.endogenous_steady_state
    exogenous = context.results.model_results[1].trends.exogenous_steady_state
    exogenous_det = context.results.model_results[1].trends.exogenous_det_steady_state
    m = context.models[1]
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = [];
        param_names = get
        for line in readlines(io)
            elem = split(line);
            if elem[1] in keys(symboltable)
                s = symboltable[elem[1]]
                if s.type == Endogenous
                    endogenous[s.orderintype] = parse(Float64, elem[2])
                elseif s.symboltype == Exogenous
                    exogenous[s.orderintype] = parse(Float64, elem[2])
                elseif s.symboltype == ExogenousDeterministic
                    exogenous_det[s.orderintype] = parse(Float64, elem[2])
                end
            end
        end
    end
    if !isempty(m.set_auxiliary_variables!)
        Base.invokelatest(m.set_auxiliary_variables!.set_auxiliary_variables!,
                          endogenous,
                          exogenous,
                          parameters)
    end
end

function isempty(m::Module)
    return names(m)[1] == :anonymous
end
