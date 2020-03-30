module Dynare

include("model.jl")
include("parser/DynareParser.jl")
include("parser/DynarePreprocessor.jl")

export parser, dynare_preprocess
end # module
