@enum PRIOR_TYPE parameter standard_deviaton correlation 

mutable struct prior
    boundaries::String
    domain
    init
    interval
    jscale
    mean
    median
    mode
    name
    name1
    name2
    parameter_type
    regimes
    shape
    shift
    stdev
    subsample
    truncate
    variance
    prior() = new()
end

function parse_prior(context, field)
    p = prior()
    p.parameter_type = parameter
    _parse(field)
end

function parse_std_prior(context, field)
    p = prior()
    p.parameter_type = parameter
    _parse(field)
end

function parse_corr_prior(context, field)
    p = prior()
    p.parameter_type = parameter
    _parse(field)
end

function _parse(field)
    for (k,v) in field
        if k == "statementName"
            nothing
        elseif k == "options"
            for (k1, v1) in v
                setfield!(p, Symbol(k1), v1)
            end
        else
            setfield!(p, Symbol(k), v)
        end
    end
    return p
end

testcase = [
Dict("statementName" => "prior", "name" => "alpha", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.356, "stdev" => 0.02)), 
Dict("statementName" => "prior", "name" => "beta", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.993, "stdev" => 0.002)), 
Dict("statementName" => "prior", "name" => "rho", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.129, "stdev" => 0.223)), 
Dict("statementName" => "prior", "name" => "delta", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.01, "stdev" => 0.005)), 
Dict("statementName" => "prior", "name" => "theta", "subsample" => "", "shape" => "normal", "options" => Dict("mean" => 3.0, "stdev" => 1.0)), 
Dict("statementName" => "prior", "name" => "tau", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.03, "stdev" => 0.01))
]

for c in testcase
    @show parse_prior([], c)
end  
