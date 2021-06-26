using SnoopCompile

#f(v::Vector{Vector{Any}}) = Vector{Vector{Int64}}(v)

function g(w)
    vv = zeros(Int64, 3)
    for (i, v) in enumerate(w)
        for (i, j) in enumerate(v)
            vv[i] = j::Int64
        end
        y = vv
    end
end


v = Vector{Any}([1, 2, 4])
w = Vector{Any}([v, v, v])

tinf = @snoopi_deep g(w)


itrigs = inference_triggers(tinf)

