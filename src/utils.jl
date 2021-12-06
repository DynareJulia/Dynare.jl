function rename_field!(d::Dict{String,Any}, name1::String, name2::String)
    d[name2] = d[name1]
    delete!(d, name1)
    return d
end

function check_isfinite(x::AbstractArray)
    if any(!isfinite, x)
        i = findall(!isfinite, x)
        throw(ErrorException("$i"))
    end
end
