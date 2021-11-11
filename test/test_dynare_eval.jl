struct Context
    symboltable::Any
    work::Any
    options::Dict
end

struct Work
    params::Any
end


function dynare_parse_eval(s::String, context::Context)
    @show s
    e = Meta.parse(s)
    @show e
    e = dynare_eval(e, context)
    return eval(e)
end

function dynare_eval(expr::Expr, context::Context)
    @show expr.head
    for (i, a) in enumerate(expr.args)
        @show a, typeof(a)
        expr.args[i] = dynare_eval(a, context)
    end
    return expr
end

function dynare_eval(s::Symbol, context)
    symboltable = context.symboltable
    ks = keys(symboltable)
    params = context.work.params
    ss = string(s)
    if ss in ks
        st = symboltable[ss]
        if st.type == Dynare.Parameter
            s = params[st.orderintype]
        end
    end
    return s
end

function dynare_eval(x::Real, context)
    return x
end

function dynare_eval(s::String, context)
    return s
end

function dynare_eval(q::QuoteNode, context)
    return q
end

function load_params!(a, b)
    @show a, b
end

context = Context(Dict(), Work([]), Dict())

dynare_parse_eval("context.options[\"LRE_solver\"] = \"cycle_reduction\"\r", context)
