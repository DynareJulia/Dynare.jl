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
