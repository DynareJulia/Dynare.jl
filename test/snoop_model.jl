using JSON
using SnoopCompile
using Test

include("../src/dynare_containers.jl")

get_model_info(field::Dict{String,Any}) = ModelInfo(
    field["lead_lag_incidence"]::Vector{Any},
    field["nstatic"]::Int64,
    field["nfwrd"]::Int64,
    field["npred"]::Int64,
    field["nboth"]::Int64,
    field["nsfwrd"]::Int64,
    field["nspred"]::Int64,
    field["ndynamic"]::Int64,
    field["maximum_endo_lag"]::Int64,
    field["maximum_endo_lead"]::Int64,
    field["maximum_exo_lag"]::Int64,
    field["maximum_exo_lead"]::Int64,
    field["maximum_exo_det_lag"]::Int64,
    field["maximum_exo_det_lead"]::Int64,
    field["maximum_lag"]::Int64,
    field["maximum_lead"]::Int64,
    field["orig_maximum_endo_lag"]::Int64,
    field["orig_maximum_endo_lead"]::Int64,
    field["orig_maximum_exo_lag"]::Int64,
    field["orig_maximum_exo_lead"]::Int64,
    field["orig_maximum_exo_det_lag"]::Int64,
    field["orig_maximum_exo_det_lead"]::Int64,
    max(
        field["orig_maximum_lag"]::Int64,
        field["orig_maximum_lag_with_diffs_expanded"]::Int64,
    ),
    field["orig_maximum_lead"]::Int64,
    field["NNZDerivatives"]::Vector{Any},
)
function load_dynare_function(modname::String, compileoption::Bool)::Module
    if compileoption
        fun = readlines(modname * ".jl")
        return (eval(Meta.parse(join(fun, "\n"))))
    else
        push!(LOAD_PATH, dirname(modname))
        name = basename(modname)
        eval(Meta.parse("using " * name))
        pop!(LOAD_PATH)
        return (eval(Symbol(name)))
    end
end

lead_lag_incidence = [
    1 0 2 0 3 0 4 0 5 0
    0 6 7 8 9 10 0 11 12 13
    14 0 15 0 16 0 17 18 0 19
]

json = """
{\"model_info\": {\"lead_lag_incidence\":
[[ 1, 0, 14], [0, 6, 0], [2, 7, 15], [0, 8, 0], [3, 9, 16], [0, 10, 0], [4, 0, 17], [0, 11, 18], [5, 12, 0], [0, 13, 19]],
\"nstatic\": 3, 
\"nfwrd\": 1, 
\"npred\": 2, 
\"nboth\": 4, 
\"nsfwrd\": 6, 
\"nspred\": 5, 
\"ndynamic\": 7, 
\"maximum_endo_lag\": 1, 
\"maximum_endo_lead\": 1, 
\"maximum_exo_lag\": 0, 
\"maximum_exo_lead\": 0, 
\"maximum_exo_det_lag\": 0, 
\"maximum_exo_det_lead\": 0, 
\"maximum_lag\": 1, 
\"maximum_lead\": 1, 
\"orig_maximum_endo_lag\": 1,
\"orig_maximum_endo_lead\": 1,
\"orig_maximum_exo_lag\": 0,
\"orig_maximum_exo_lead\": 0,
\"orig_maximum_exo_det_lag\": 0,
\"orig_maximum_exo_det_lead\": 0,
\"orig_maximum_lag\": 1,
\"orig_maximum_lead\": 1,
\"orig_maximum_lag_with_diffs_expanded\": 1,
\"NNZDerivatives\": [26, -1, -1]}}
"""

modfile = JSON.parse(json)

@show modfile["model_info"]

model_info = get_model_info(modfile["model_info"])

NNZDerivatives = [n for n::Int64 in model_info.NNZDerivatives]

modfilename = "models/example1/example1"
endo_nbr = 10
exo_nbr = 1
exo_det_nbr = 0
param_nbr = 1

function f(
    modfilename,
    endo_nbr,
    lead_lag_incidence,
    exo_nbr,
    lagged_exo_nbr,
    exo_det_nbr,
    param_nbr,
)

    model = Model(
        modfilename,
        endo_nbr,
        lead_lag_incidence,
        exo_nbr,
        lagged_exo_nbr,
        exo_det_nbr,
        param_nbr,
        model_info.maximum_endo_lag,
        model_info.maximum_endo_lead,
        model_info.maximum_exo_lag,
        model_info.maximum_exo_lead,
        model_info.maximum_exo_det_lag,
        model_info.maximum_exo_det_lead,
        model_info.maximum_lag,
        model_info.maximum_lead,
        model_info.orig_maximum_endo_lag,
        model_info.orig_maximum_endo_lead,
        model_info.orig_maximum_exo_lag,
        model_info.orig_maximum_exo_lead,
        model_info.orig_maximum_exo_det_lag,
        model_info.orig_maximum_exo_det_lead,
        model_info.orig_maximum_lag,
        model_info.orig_maximum_lead,
        NNZDerivatives,
        true,
    )
end

tinf = @snoopi_deep f(
    modfilename,
    endo_nbr,
    model_info.lead_lag_incidence,
    exo_nbr,
    0,
    exo_det_nbr,
    param_nbr,
)

itrigs = inference_triggers(tinf)
mtrigs = accumulate_by_source(Method, itrigs)
