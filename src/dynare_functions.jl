module DFunctions

using RuntimeGeneratedFunctions
using StatsFuns
using TimeDataFrames

RuntimeGeneratedFunctions.init(@__MODULE__)

function load_set_dynamic_auxiliary_variables(modelname::String)
    source = []
    functionstart = false
    for line in readlines("$(modelname)DynamicSetAuxiliarySeries.jl", keep = true)
        if startswith(line, "function")
            functionstart = true
        end
        if functionstart
            push!(source, line)
            if startswith(line, "end")
#                push!(source, "end")
                break
            end
        end
    end
    exp1 = Meta.parse(join(source, "\n"))
    #    convert_expression(exp1)
    return (@RuntimeGeneratedFunction(exp1))
end


#=
is_ds_var(e::Expr) =  e.head == :. && e.args[1] == :ds
is_ds_arg(e::Expr) = e.head == :call && typeof(e.args[1]) == Expr && is_ds_var(e.args[1])

function convert_expression(e)
    if is_ds_arg(e)
        k = e.args[2]
        e.args[2] = e.args[1]
        e.args[1] = (k < 0) ? :lag : :lead 
        if abs(k) > 1
            push!(e.args, abs(k))
        end
    end
    for a in e.args
        if typeof(a) == Expr
            convert_expression(a)
        end
    end
end
=#
end # end module
