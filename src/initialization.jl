function histval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    histval = zeros(m.orig_maximum_lag, m.endogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]::String].orderintype::Int64
        l = m.orig_maximum_lag - v["lag"]::Int64
        histval[l, k] = dynare_parse_eval(v["value"]::String, context)
    end
    context.work.histval .= histval
end

function param_init!(context::Context, field::Dict{String, Any})
    params = context.work.params
    symboltable::Dict{String, DynareSymbol} = context.symboltable
    s = symboltable[field["name"]::String]
    k = s.orderintype::Int64
    params[k] = Float64(dynare_parse_eval(field["value"]::String, context))
end

function initval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    endo_steady_state = zeros(m.endogenous_nbr)
    exo_steady_state = zeros(m.exogenous_nbr)
    exo_det_steady_state = zeros(m.exogenous_deterministic_nbr)
    for v in field["vals"]
        s = symboltable[v["name"]::String]
        typ = s.symboltype
        k = s.orderintype::Int64
        value = dynare_parse_eval(v["value"]::String, context)
        if typ == :Endogenous
            endo_steady_state[k] = value
        elseif typ == :Exogenous
            exo_steady_state[k] = value
        elseif typ == :ExogenousDeterministic
            exo_det_steady_state[k] = value
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
        k =  symboltable[v["name"]::String].orderintype::Int64
        Sigma[k, k] = dynare_parse_eval(v["variance"]::String, context)::Float64
    end
end

function set_stderr!(Sigma::Matrix{Float64}, stderr::Vector{Any}, symboltable::SymbolTable)
    for s in stderr
        k =  symboltable[s["name"]::String].orderintype::Int64
        x = dynare_parse_eval(s["stderr"]::String, context)::Float64
        Sigma[k, k] = x*x
    end
end

function set_covariance!(Sigma::Matrix{Float64}, covariance::Vector{Any}, symboltable::SymbolTable)
    for c in covariance
        k1 =  symboltable[c["name"]::String].orderintype::Int64
        k2 =  symboltable[c["name2"]::String].orderintype::Int64
        Sigma[k1, k2] = dynare_parse_eval(c["covariance"]::String, context)::Float64
        Sigma[k2, k1] = Sigma[k1, k2]
    end
end

function set_correlation!(Sigma::Matrix{Float64}, correlation::Vector{Any}, symboltable::SymbolTable)
    for c in correlation
        k1 =  symboltable[c["name"]::String].orderintype::Int64
        k2 =  symboltable[c["name2"]::String].orderintype::Int64
        corr = dynare_parse_eval(c["correlation"]::String, context)::Float64
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
                    k = s.orderintype::Int64
                    parameters[k] = parse(Float64, elem[2])
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
                k = s.orderintype::Int64
                if s.symboltype == Endogenous
                    endogenous[k] = parse(Float64, elem[2])
                elseif s.symboltype == Exogenous
                    exogenous[k] = parse(Float64, elem[2])
                elseif s.symboltype == ExogenousDeterministic
                    exogenous_det[k] = parse(Float64, elem[2])
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
