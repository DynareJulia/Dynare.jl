macro limits(block)
    for a in block.args
        if typeof(a) == Expr
            t = eval(a)
            limits!(Symbol(t[1]), min = t[2], max = t[3])
        end
    end
    context.work.limits = limits
end

function limits!(s; context=context, min=-Inf, max=Inf)
    context.work.limits[Symbol(s)] = (max=dynare_eval(max, context, [], []), min=dynare_eval(min, context, [], []))
end
