function histval!(context, field)
    symboltable = context.symboltable
    m = context.models[1]
    histval = zeros(m.orig_maximum_lag, m.endogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]].orderintype
        l = m.orig_maximum_lag - v["lag"]
        histval[l, k] = Meta.parse(v["value"])
    end
    context.work.histval = histval
end

function param_init!(context, field)
    params = context.work.params;
    symboltable = context.symboltable
    s = symboltable[field["name"]]
    k = s.orderintype
    params[k] = eval(Meta.parse(field["value"]))
end

function initval!(context, field)
    symboltable = context.symboltable
    m = context.models[1]
    endo_steady_state = zeros(m.endogenous_nbr)
    exo_steady_state = zeros(m.exogenous_nbr)
    exo_det_steady_state = zeros(m.exogenous_deterministic_nbr)
    for v in field["vals"]
        s = symboltable[v["name"]]
        typ = s.type
        k = s.orderintype
        if typ == :Endogenous
            endo_steady_state[k] = Meta.parse(v["value"])
        elseif typ == :Exogenous
            exo_steady_state[k] = Meta.parse(v["value"])
        elseif typ == :ExogenousDeterministic
            exo_det_steady_state[k] = Meta.parse(v["value"])
        else
            throw(error("$(v["name"]) can't be set in INITVAL"))
        end
    end
    params = context.work.params
    Base.invokelatest(
        m.set_auxiliary_variables!.set_auxiliary_variables!(
            endogenous_steady_state,
            exogenous_steady_state,
            params)
    )
    trends = context.results.results_model[1].trends 
    trends.endogenous_steady_state = endo_steady_state
    trends.exogenous_steady_state = exo_steady_state
    trends.exogenous_det_steady_state = exo_det_steady_state
end

function shocks!(context, field)
    Sigma = context.models[1].Sigma_e
    symboltable = context.symboltable
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
        x = eval(Meta.parse(s["stderr"]))
        Sigma[k, k] = x*x
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
                if s.type == Parameter
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
                elseif s.type == Exogenous
                    exogenous[s.orderintype] = parse(Float64, elem[2])
                elseif s.type == ExogenousDeterministic
                    exogenous_det[s.orderintype] = parse(Float64, elem[2])
                end
            end
        end
    end
    Base.invokelatest(m.set_auxiliary_variables!,
                      endogenous,
                      exogenous,
                      parameters)
end
