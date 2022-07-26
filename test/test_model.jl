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

model = Model(
    modfilename,
    endo_nbr,
    model_info.lead_lag_incidence,
    exo_nbr,
    0,
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
    #              commandlineoptions.compilemodule
    true,
)

@test model.endogenous_nbr == endo_nbr
@test model.exogenous_nbr == exo_nbr
@test model.lagged_exogenous_nbr == 0
@test model.exogenous_deterministic_nbr == exo_det_nbr
@test model.parameter_nbr == param_nbr
@test model.lead_lag_incidence == lead_lag_incidence
@test model.n_static == 3
@test model.n_fwrd == 2
@test model.n_bkwrd == 1
@test model.n_both == 4
@test model.n_states == 5
@test model.DErows1 == 1:7
@test model.DErows2 == 8:11
@test model.n_dyn == 11
@test model.i_static == [2, 4, 6]
@test model.i_dyn == [1, 3, 5, 7, 8, 9, 10]
@test model.i_bkwrd == [9]
@test model.i_bkwrd_b == [1, 3, 5, 7, 9]
@test model.i_bkwrd_ns == [1, 2, 3, 4, 6]
@test model.i_fwrd == [8, 10]
@test model.i_fwrd_b == [1, 3, 5, 7, 8, 10]
@test model.i_fwrd_ns == [1, 2, 3, 4, 5, 7]
@test model.i_both == [1, 3, 5, 7]
@test model.i_non_states == [2, 4, 6, 8, 10]
@test model.p_static == [6, 8, 10]
@test model.p_bkwrd == [5]
@test model.p_bkwrd_b == [1, 2, 3, 4, 5]
@test model.p_fwrd == [18, 19]
@test model.p_fwrd_b == [14, 15, 16, 17, 18, 19]
@test model.p_both_b == [1, 2, 3, 4]
@test model.p_both_f == [14, 15, 16, 17]
@test model.i_current == [2, 3, 4, 5, 6, 8, 9, 10]
@test model.p_current == 6:13
@test model.n_current == 8
@test model.i_current_ns == [2, 3, 5, 6, 7]
@test model.p_current_ns == [7, 9, 11, 12, 13]
@test model.n_current_ns == 5
@test model.icolsD == [5, 6, 7, 8, 9, 10, 11]
@test model.jcolsD == [12, 14, 15, 16, 17, 18, 19]
@test model.icolsE == [1, 2, 3, 4, 5, 7, 8, 10, 11]
@test model.jcolsE == [1, 2, 3, 4, 5, 7, 9, 11, 13]
@test model.colsUD == [1, 2, 3, 4]
@test model.colsUE == [6, 7, 8, 9]
@test model.i_cur_fwrd == [2, 3, 5, 6]
@test model.n_cur_fwrd == 4
@test model.p_cur_fwrd == [7, 9, 11, 13]
@test model.i_cur_bkwrd == [5]
@test model.n_cur_bkwrd == 1
@test model.p_cur_bkwrd == [12]
@test model.i_cur_both == [2, 3]
@test model.n_cur_both == 2
@test model.p_cur_both == [7, 9]
@test model.gx_rows == 6:11
@test model.hx_rows == 1:5
@test model.i_current_exogenous == [20]
@test model.i_lagged_exogenous == []
@test model.serially_correlated_exogenous == []
#@test model.Sigma_e == 
@test model.maximum_endo_lag == 1
@test model.maximum_endo_lead == 1
@test model.maximum_exo_lag == 0
@test model.maximum_exo_lead == 0
@test model.maximum_exo_det_lag == 0
@test model.maximum_exo_det_lead == 0
@test model.maximum_lag == 1
@test model.maximum_lead == 1
@test model.orig_maximum_endo_lag == 1
@test model.orig_maximum_endo_lead == 1
@test model.orig_maximum_exo_lag == 0
@test model.orig_maximum_exo_lead == 0
@test model.orig_maximum_exo_det_lag == 0
@test model.orig_maximum_exo_det_lead == 0
@test model.orig_maximum_lag == 1
@test model.orig_maximum_lead == 1
@test model.dynamic_indices == [1, 3, 5, 7, 8, 9, 10]
@test model.current_dynamic_indices == [3, 5, 8, 9, 10]
@test model.forward_indices_d == [5, 7]
@test model.backward_indices_d == [6]
@test model.current_dynamic_indices_d == [2, 3, 5, 6, 7]
@test model.exogenous_indices == [20]
@test model.NNZDerivatives == [26, -1, -1]
