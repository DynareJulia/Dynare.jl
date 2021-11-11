using JSON
using SnoopCompile
using Test

include("../src/dynare_containers.jl")

function f(k::String, s::SymbolTable)
    s1 = s[k]
end

struct DynareSymbol
    longname::String
    texname::String
    symboltype::SymbolType
    orderintype::Integer
end

v1 = DynareSymbol("longname1", "texname1", Endogenous, 1)
v2 = DynareSymbol("longname2", "texname2", Endogenous, 2)

st = Dict("v1" => v1, "v2" => v2)

tinf = @snoopi_deep f("v2", st)
itrigs = inference_triggers(tinf)
mtrigs = accumulate_by_source(Method, itrigs)
