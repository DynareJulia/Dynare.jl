function rename_field!(d::Dict{String, Any}, name1::String, name2::String)
    d[name2] = d[name1]
    delete!(d, name1)
    return d
end
