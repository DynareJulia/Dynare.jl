import Base
using Distributions
using LinearRationalExpectations
using RuntimeGeneratedFunctions
using Suppressor
using StatsFuns
using AxisArrayTables
using KalmanFilterTools: KalmanLikelihoodWs

export Context,
    DynareSymbol, Model, ModelResults, Results, Simulation, SymbolType, Work, Trends

RuntimeGeneratedFunctions.init(@__MODULE__)

@enum SymbolType Endogenous Exogenous ExogenousDeterministic Parameter DynareFunction
@enum EstimatedParameterType EstParameter EstSDShock EstSDMeasurement EstVarShock EstVarMeasurement EstCorrShock EstCorrMeasurement

mutable struct ModFileInfo
    endval_is_reset::Bool
    has_auxiliary_variables::Bool
    has_calib_smoother::Bool
    has_check::Bool
    has_deterministic_trend::Bool
    has_dynamic_file::Bool
    has_endval::Bool
    has_histval::Bool
    has_histval_file::Bool
    has_initval::Bool
    has_initval_file::Bool
    has_planner_objective::Bool
    has_perfect_foresight_setup::Bool
    has_perfect_foresight_solver::Bool
    has_ramsey_model::Bool
    has_shocks::Bool
    has_static_file::Bool
    has_steadystate_file::Bool
    has_stoch_simul::Bool
    has_trends::Bool
    initval_is_reset::Bool
    modfilepath::String
   
    function ModFileInfo(modfilepath_arg::String)
        endval_is_reset = false
        has_auxiliary_variables = false
        has_calib_smoother = false
        has_check = false
        has_deterministic_trend = false
        has_dynamic_file = false
        has_endval = false
        has_histval = false
        has_histval_file = false
        has_initval = false
        has_initval_file = false
        has_planner_objective = false
        has_perfect_foresight_setup = false
        has_perfect_foresight_solver = false
        has_ramsey_model = false
        has_shocks = false
        has_static_file = false
        has_steadystate_file = false
        has_stoch_simul = false
        has_trends = false
        initval_is_reset = false
        modfilepath = modfilepath_arg
        new(
            endval_is_reset,
            has_auxiliary_variables,
            has_calib_smoother,
            has_check,
            has_deterministic_trend,
            has_dynamic_file,
            has_endval,
            has_histval,
            has_histval_file,
            has_initval,
            has_initval_file,
            has_planner_objective,
            has_perfect_foresight_setup,
            has_perfect_foresight_solver,
            has_ramsey_model,
            has_shocks,
            has_static_file,
            has_steadystate_file,
            has_stoch_simul,
            has_trends,
            initval_is_reset,
            modfilepath,
        )
    end
end

Base.show(io::IO, mf::ModFileInfo) = show_field_value(mf)

struct Model
    endogenous_nbr::Int
    exogenous_nbr::Int
    lagged_exogenous_nbr::Int
    exogenous_deterministic_nbr::Int
    parameter_nbr::Int
    original_endogenous_nbr::Int
    lead_lag_incidence::Matrix{Int}
    n_static::Int
    n_fwrd::Int
    n_bkwrd::Int
    n_both::Int
    n_states::Int
    DErows1::Vector{Int}
    DErows2::Vector{Int}
    n_dyn::Int
    i_static::Vector{Int}
    i_dyn::Array{Int,1}
    i_bkwrd::Vector{Int}
    i_bkwrd_b::Vector{Int}
    i_bkwrd_ns::Vector{Int}
    i_fwrd::Vector{Int}
    i_fwrd_b::Vector{Int}
    i_fwrd_ns::Vector{Int}
    i_both::Vector{Int}
    i_non_states::Vector{Int}
    p_static::Vector{Int}
    p_bkwrd::Vector{Int}
    p_bkwrd_b::Vector{Int}
    p_fwrd::Vector{Int}
    p_fwrd_b::Vector{Int}
    p_both_b::Vector{Int}
    p_both_f::Vector{Int}
    i_current::Vector{Int}
    p_current::Vector{Int}
    n_current::Int
    i_current_ns::Vector{Int}
    p_current_ns::Vector{Int}
    n_current_ns::Int
    icolsD::Vector{Int}
    jcolsD::Vector{Int}
    icolsE::Vector{Int}
    jcolsE::Vector{Int}
    colsUD::Vector{Int}
    colsUE::Vector{Int}
    i_cur_fwrd::Vector{Int}
    n_cur_fwrd::Int
    p_cur_fwrd::Vector{Int}
    i_cur_bkwrd::Vector{Int}
    n_cur_bkwrd::Int
    p_cur_bkwrd::Vector{Int}
    i_cur_both::Vector{Int}
    n_cur_both::Int
    p_cur_both::Vector{Int}
    gx_rows::Vector{Int}
    hx_rows::Vector{Int}
    i_current_exogenous::Vector{Int}
    i_lagged_exogenous::Vector{Int}
    serially_correlated_exogenous::Vector{Int}
    Sigma_e::Matrix{Float64}
    maximum_endo_lag::Int
    maximum_endo_lead::Int
    maximum_exo_lag::Int
    maximum_exo_lead::Int
    maximum_exo_det_lag::Int
    maximum_exo_det_lead::Int
    maximum_lag::Int
    maximum_lead::Int
    orig_maximum_endo_lag::Int
    orig_maximum_endo_lead::Int
    orig_maximum_exo_lag::Int
    orig_maximum_exo_lead::Int
    orig_maximum_exo_det_lag::Int
    orig_maximum_exo_det_lead::Int
    orig_maximum_lag::Int
    orig_maximum_lead::Int
    dynamic_indices::Vector{Int}
    current_dynamic_indices::Vector{Int}
    forward_indices_d::Vector{Int}
    backward_indices_d::Vector{Int}
    current_dynamic_indices_d::Vector{Int}
    exogenous_indices::Vector{Int}
    NNZDerivatives::Vector{Int}
    auxiliary_variables::Vector{Dict{String,Any}}
    mcps::Vector{Tuple{Int,Int,String,String}}
    dynamic_g1_sparse_rowval::Vector{Int}
    dynamic_g1_sparse_colptr::Vector{Int}
    dynamic_g2_sparse_indices::Vector{Vector{Int}}
    static_g1_sparse_rowval::Vector{Int}
    static_g1_sparse_colptr::Vector{Int}    
    dynamic_tmp_nbr::Vector{Int}
    static_tmp_nbr::Vector{Int}
    ids::LRE.Indices
end


# purely backward model
function assemble_lead_lag_incidence_1!(
    lead_lag_incidence_matrix::Matrix{Int},
    backward_indices_d::Vector{Int},
    current_dynamic_indices::Vector{Int},
    i_backward_in_current::Vector{Int},
    i_bkwrd::Vector{Int},
    i_bkwrd_ns::Vector{Int},
    i_cur_bkwrd::Vector{Int},
    i_current::Vector{Int},
    i_current_ns::Vector{Int},
    i_dyn::Vector{Int},
    i_non_states::Vector{Int},
    i_static::Vector{Int},
    p_bkwrd::Vector{Int},
    p_cur_bkwrd::Vector{Int},
    p_current::Vector{Int},
    p_current_ns::Vector{Int},
    p_static::Vector{Int},
    lead_lag_incidence::Vector{Vector{Int}},
    endogenous_nbr::Int,
)
    max_lead_lag_incidence = 0
    n_dyn1 = 0
    n_current = 0
    n_bkwrd = 0
    n_static = 0

    for (i, v) in enumerate(lead_lag_incidence)
        for j = 1:2
            lead_lag_incidence_matrix[j, i] = v[j]
        end
        if v[1] > 0
            n_dyn1 += 1
            push!(i_dyn, i)
            push!(i_bkwrd_ns, n_dyn1)
            max_lead_lag_incidence = max(v[1], max_lead_lag_incidence)
            n_bkwrd += 1
            push!(i_bkwrd, i)
            push!(p_bkwrd, v[1])
            push!(backward_indices_d, n_dyn1)
        else
            push!(i_non_states, i)
            n_static += 1
            push!(i_static, i)
            push!(p_static, v[2])
        end
        if v[2] > 0
            n_current += 1
            push!(i_current, i)
            push!(p_current, v[2])
            max_lead_lag_incidence = max(v[2], max_lead_lag_incidence)
            if v[1] > 0
                push!(i_current_ns, n_dyn1)
                push!(p_current_ns, v[2])
                push!(current_dynamic_indices, i)
                push!(i_cur_bkwrd, n_bkwrd)
                push!(p_cur_bkwrd, v[2])
                push!(i_backward_in_current, n_bkwrd)
            end
        end
    end
    return (max_lead_lag_incidence, n_bkwrd, n_current, n_static)

end

# backward and forward model
function assemble_lead_lag_incidence_3!(
    lead_lag_incidence_matrix::Matrix{Int},
    backward_indices_d::Vector{Int},
    current_dynamic_indices::Vector{Int},
    i_backward_in_current::Vector{Int},
    i_bkwrd::Vector{Int},
    i_bkwrd_b::Vector{Int},
    i_bkwrd_ns::Vector{Int},
    i_both::Vector{Int},
    i_both_b::Vector{Int},
    i_both_f::Vector{Int},
    i_cur_bkwrd::Vector{Int},
    i_cur_both::Vector{Int},
    i_cur_fwrd_b::Vector{Int},
    i_current::Vector{Int},
    i_current_ns::Vector{Int},
    i_dyn::Vector{Int},
    i_forward_in_current::Vector{Int},
    i_fwrd::Vector{Int},
    i_fwrd_b::Vector{Int},
    i_fwrd_ns::Vector{Int},
    i_non_states::Vector{Int},
    i_static::Vector{Int},
    forward_indices_d::Vector{Int},
    p_bkwrd::Vector{Int},
    p_bkwrd_b::Vector{Int},
    p_both_b::Vector{Int},
    p_both_f::Vector{Int},
    p_cur_bkwrd::Vector{Int},
    p_current::Vector{Int},
    p_current_ns::Vector{Int},
    p_fwrd::Vector{Int},
    p_fwrd_b::Vector{Int},
    p_static::Vector{Int},
    p_cur_fwrd_b::Vector{Int},
    p_cur_both::Vector{Int},
    lead_lag_incidence::Vector{Vector{Int}},
    endogenous_nbr::Int,
)
    max_lead_lag_incidence = 0
    n_dyn1 = 0
    n_dyn2 = 0
    n_bkwrd_b = 0
    n_current = 0
    n_fwrd_b = 0
    n_both = 0
    n_bkwrd = 0
    n_fwrd = 0
    n_static = 0

    for (i, v) in enumerate(lead_lag_incidence)
        for j = 1:3
            lead_lag_incidence_matrix[j, i] = v[j]
        end
        if v[1] > 0
            n_dyn1 += 1
            n_bkwrd_b += 1
            push!(i_bkwrd_b, i)
            push!(p_bkwrd_b, v[1])
            push!(i_dyn, i)
            push!(i_bkwrd_ns, n_dyn1)
            max_lead_lag_incidence = max(max(v[1], v[3]), max_lead_lag_incidence)
            if v[3] > 0
                n_dyn2 += 1
                n_fwrd_b += 1
                n_both += 1
                push!(i_both, i)
                push!(i_both_b, n_bkwrd_b)
                push!(i_both_f, n_fwrd_b)
                push!(i_fwrd_b, i)
                push!(i_fwrd_ns, n_dyn1)
                push!(p_both_b, v[1])
                push!(p_both_f, v[3])
                push!(p_fwrd_b, v[3])
            else
                n_bkwrd += 1
                push!(i_bkwrd, i)
                push!(p_bkwrd, v[1])
                push!(backward_indices_d, n_dyn1)
            end
        else
            push!(i_non_states, i)
            if v[3] > 0
                n_dyn1 += 1
                n_dyn2 += 1
                n_fwrd += 1
                n_fwrd_b += 1
                push!(i_dyn, i)
                push!(i_fwrd, i)
                push!(p_fwrd, v[3])
                push!(i_fwrd_b, i)
                push!(p_fwrd_b, v[3])
                push!(i_fwrd_ns, n_dyn1)
                push!(forward_indices_d, n_dyn1)
                max_lead_lag_incidence = max(v[3], max_lead_lag_incidence)
            else
                n_static += 1
                push!(i_static, i)
                push!(p_static, v[2])
            end
        end
        if v[2] > 0
            n_current += 1
            push!(i_current, i)
            push!(p_current, v[2])
            max_lead_lag_incidence = max(v[2], max_lead_lag_incidence)
            if v[1] > 0 || v[3] > 0
                push!(i_current_ns, n_dyn1)
                push!(p_current_ns, v[2])
                push!(current_dynamic_indices, i)
            end
            if v[1] > 0
                if v[3] == 0
                    push!(i_cur_bkwrd, n_bkwrd_b)
                    push!(p_cur_bkwrd, v[2])
                end
                push!(i_backward_in_current, n_bkwrd_b)
            end
            if v[3] > 0
                push!(i_cur_fwrd_b, n_fwrd_b)
                push!(p_cur_fwrd_b, v[2])
                push!(i_forward_in_current, n_fwrd_b)
            end
            if v[1] > 0 && v[3] > 0
                push!(i_cur_both, n_both)
                push!(p_cur_both, v[2])
            end

        end
    end
    return (max_lead_lag_incidence, n_bkwrd, n_both, n_current, n_fwrd, n_static)

end

function Model(
    modfilename::String,
    modfileinfo::ModFileInfo,
    endogenous_nbr::Int,
    lead_lag_incidence::Vector{Vector{Int}},
    exogenous_nbr::Int,
    lagged_exogenous_nbr::Int,
    exogenous_deterministic_nbr::Int,
    parameter_nbr::Int,
    orig_endo_nbr::Int,
    aux_vars::Vector{Any},
    maximum_endo_lag::Int,
    maximum_endo_lead::Int,
    maximum_exo_lag::Int,
    maximum_exo_lead::Int,
    maximum_exo_det_lag::Int,
    maximum_exo_det_lead::Int,
    maximum_lag::Int,
    maximum_lead::Int,
    orig_maximum_endo_lag::Int,
    orig_maximum_endo_lead::Int,
    orig_maximum_exo_lag::Int,
    orig_maximum_exo_lead::Int,
    orig_maximum_exo_det_lag::Int,
    orig_maximum_exo_det_lead::Int,
    orig_maximum_lag::Int,
    orig_maximum_lead::Int,
    NNZDerivatives::Vector{Int},
    compileoption::Bool,
    dynamic_g1_sparse_rowval::Vector{Int},
    dynamic_g1_sparse_colptr::Vector{Int},
    dynamic_g2_sparse_indices::Vector{Vector{Int}},
    static_g1_sparse_rowval::Vector{Int},
    static_g1_sparse_colptr::Vector{Int},
    dynamic_tmp_nbr::Vector{Int},
    static_tmp_nbr::Vector{Int}
)
    i_static = Vector{Int}(undef, 0)
    p_static = similar(i_static)
    i_dyn = similar(i_static)
    i_bkwrd = similar(i_static)
    i_bkwrd_b = similar(i_static)
    i_bkwrd_ns = similar(i_static)
    p_bkwrd = similar(i_static)
    p_bkwrd_b = similar(i_static)
    i_fwrd = similar(i_static)
    i_fwrd_b = similar(i_static)
    i_fwrd_ns = similar(i_static)
    p_fwrd = similar(i_static)
    p_fwrd_b = similar(i_static)
    i_both = similar(i_static)
    i_non_states = similar(i_static)
    i_both_b = similar(i_static)
    p_both_b = similar(i_static)
    i_both_f = similar(i_static)
    p_both_f = similar(i_static)
    i_current = similar(i_static)
    p_current = similar(i_static)
    i_current_ns = similar(i_static)
    p_current_ns = similar(i_static)
    i_cur_fwrd_b = similar(i_static)
    p_cur_fwrd_b = similar(i_static)
    i_cur_bkwrd = similar(i_static)
    p_cur_bkwrd = similar(i_static)
    i_cur_both = similar(i_static)
    p_cur_both = similar(i_static)
    i_backward_in_current = similar(i_static)
    i_forward_in_current = similar(i_static)
    icolsD = similar(i_static)
    jcolsD = similar(i_static)
    icolsE = similar(i_static)
    jcolsE = similar(i_static)
    colsUD = similar(i_static)
    colsUE = similar(i_static)
    n_dyn = similar(i_static)
    DErows1 = similar(i_static)
    DErows2 = similar(i_static)
    gx_rows = similar(i_static)
    hx_rows = similar(i_static)
    current_dynamic_indices = similar(i_static)
    forward_indices_d = similar(i_static)
    backward_indices_d = similar(i_static)
    current_dynamic_indices_d = similar(i_static)
    if maximum_endo_lag == 1
        if maximum_endo_lead == 1
            # backward and forward
            lead_lag_incidence_matrix = zeros(Int, 3, endogenous_nbr)
            (max_lead_lag_incidence, n_bkwrd, n_both, n_current, n_fwrd, n_static) =
                assemble_lead_lag_incidence_3!(
                    lead_lag_incidence_matrix,
                    backward_indices_d,
                    current_dynamic_indices,
                    i_backward_in_current,
                    i_bkwrd,
                    i_bkwrd_b,
                    i_bkwrd_ns,
                    i_both,
                    i_both_b,
                    i_both_f,
                    i_cur_bkwrd,
                    i_cur_both,
                    i_cur_fwrd_b,
                    i_current,
                    i_current_ns,
                    i_dyn,
                    i_forward_in_current,
                    i_fwrd,
                    i_fwrd_b,
                    i_fwrd_ns,
                    i_non_states,
                    i_static,
                    forward_indices_d,
                    p_bkwrd,
                    p_bkwrd_b,
                    p_both_b,
                    p_both_f,
                    p_cur_bkwrd,
                    p_current,
                    p_current_ns,
                    p_fwrd,
                    p_fwrd_b,
                    p_static,
                    p_cur_fwrd_b,
                    p_cur_both,
                    lead_lag_incidence,
                    endogenous_nbr,
                )
        else
            # purely backward
            lead_lag_incidence_matrix = zeros(Int, 2, endogenous_nbr)
            (max_lead_lag_incidence, n_bkwrd, n_current, n_static) =
                assemble_lead_lag_incidence_1!(
                    lead_lag_incidence_matrix::Matrix{Int},
                    backward_indices_d::Vector{Int},
                    current_dynamic_indices::Vector{Int},
                    i_backward_in_current::Vector{Int},
                    i_bkwrd::Vector{Int},
                    i_bkwrd_ns::Vector{Int},
                    i_cur_bkwrd::Vector{Int},
                    i_current::Vector{Int},
                    i_current_ns::Vector{Int},
                    i_dyn::Vector{Int},
                    i_non_states::Vector{Int},
                    i_static::Vector{Int},
                    p_bkwrd::Vector{Int},
                    p_cur_bkwrd::Vector{Int},
                    p_current::Vector{Int},
                    p_current_ns::Vector{Int},
                    p_static::Vector{Int},
                    lead_lag_incidence::Vector{Vector{Int}},
                    endogenous_nbr::Int,
                )
            i_bkwrd_b = i_bkwrd
            p_bkwrd_b = p_bkwrd
            i_fwrd_b = i_fwrd
            p_fwrd_b = p_fwrd
        end
    elseif maximum_endo_lead == 1
        # purely forward
        assemble_lead_lag_incidence_2()
    else
        # static
        assemble_lead_lag_incidence_0()
    end
    n_bkwrd = length(i_bkwrd)
    n_fwrd = length(i_fwrd)
    n_both = length(i_both)
    n_states = n_bkwrd + n_both
    n_current = length(i_current)
    n_current_ns = length(i_current_ns)
    n_cur_fwrd_b = length(i_cur_fwrd_b)
    n_cur_bkwrd = length(i_cur_bkwrd)
    n_cur_both = length(i_cur_both)

    icolsD = vcat(i_cur_bkwrd, n_bkwrd + n_both .+ collect(1:(n_fwrd+n_both)))
    jcolsD = vcat(p_cur_bkwrd, p_fwrd_b)
    # derivatives of current values of variables that are both
    # forward and backward are included in the E matrix
    icolsE = vcat(collect(1:(n_bkwrd+n_both)), n_bkwrd + n_both .+ i_cur_fwrd_b)
    jcolsE = vcat(p_bkwrd_b, p_cur_fwrd_b)
    colsUD = i_both_b
    colsUE = n_bkwrd + n_both .+ i_both_f
    n_dyn = endogenous_nbr - n_static + n_both
    DErows1 = collect(1:(n_dyn-n_both))
    DErows2 = (n_dyn - n_both) .+ collect(1:n_both)
    gx_rows = n_bkwrd + n_both .+ collect(1:(n_fwrd+n_both))
    hx_rows = collect(1:(n_bkwrd+n_both))
    i_current_exogenous = max_lead_lag_incidence .+ (1:exogenous_nbr)
    i_lagged_exogenous = []
    Sigma_e = zeros(exogenous_nbr, exogenous_nbr)
    serially_correlated_exogenous = []
    static_indices = i_static
    current_indices = i_current
    forward_indices = i_fwrd
    both_indices = i_both
    backward_indices = i_bkwrd
    backward_number = n_bkwrd
    forward_number = n_fwrd
    both_number = n_both
    current_number = n_current
    # purely_forward_indices = [i for i in forward_indices if !(i in both_indices)]
    #    forward_indices_d = findall(in(forward_indices), i_dyn)
    #    backward_indices_d = findall(in(backward_indices), i_dyn)
    current_dynamic_indices_d = i_current_ns
    exogenous_indices = (
        backward_number + current_number + forward_number + 2 * both_number .+
        collect(1:exogenous_nbr)
    )
    mcps = Tuple{Int,String,String,String}[]
    ids = LRE.Indices(exogenous_nbr, i_fwrd_b, i_current, i_bkwrd_b, i_static)
    Model(
        endogenous_nbr,
        exogenous_nbr,
        lagged_exogenous_nbr,
        exogenous_deterministic_nbr,
        parameter_nbr,
        orig_endo_nbr,
        lead_lag_incidence_matrix,
        n_static,
        n_fwrd,
        n_bkwrd,
        n_both,
        n_states,
        DErows1,
        DErows2,
        n_dyn,
        i_static,
        i_dyn,
        i_bkwrd,
        i_bkwrd_b,
        i_bkwrd_ns,
        i_fwrd,
        i_fwrd_b,
        i_fwrd_ns,
        i_both,
        i_non_states,
        p_static,
        p_bkwrd,
        p_bkwrd_b,
        p_fwrd,
        p_fwrd_b,
        p_both_b,
        p_both_f,
        i_current,
        p_current,
        n_current,
        i_current_ns,
        p_current_ns,
        n_current_ns,
        icolsD,
        jcolsD,
        icolsE,
        jcolsE,
        colsUD,
        colsUE,
        i_cur_fwrd_b,
        n_cur_fwrd_b,
        p_cur_fwrd_b,
        i_cur_bkwrd,
        n_cur_bkwrd,
        p_cur_bkwrd,
        i_cur_both,
        n_cur_both,
        p_cur_both,
        gx_rows,
        hx_rows,
        i_current_exogenous,
        i_lagged_exogenous,
        serially_correlated_exogenous,
        Sigma_e,
        maximum_endo_lag,
        maximum_endo_lead,
        maximum_exo_lag,
        maximum_exo_lead,
        maximum_exo_det_lag,
        maximum_exo_det_lead,
        maximum_lag,
        maximum_lead,
        orig_maximum_endo_lag,
        orig_maximum_endo_lead,
        orig_maximum_exo_lag,
        orig_maximum_exo_lead,
        orig_maximum_exo_det_lag,
        orig_maximum_exo_det_lead,
        orig_maximum_lag,
        orig_maximum_lead,
        i_dyn,
        current_dynamic_indices,
        forward_indices_d,
        backward_indices_d,
        current_dynamic_indices_d,
        exogenous_indices,
        NNZDerivatives,
        aux_vars,
        mcps,
        dynamic_g1_sparse_rowval,
        dynamic_g1_sparse_colptr,
        dynamic_g2_sparse_indices,
        static_g1_sparse_rowval,
        static_g1_sparse_colptr,
        dynamic_tmp_nbr,
        static_tmp_nbr,
        ids
    )

end

Base.show(io::IO, m::Model) = show_field_value(m)

struct Simulation
    firstperiod::PeriodsSinceEpoch
    lastperiod::PeriodsSinceEpoch
    name::String
    statement::String
    data::AxisArrayTable
end

Base.show(io::IO, s::Simulation) = show_field_value(s)

mutable struct Trends
    endogenous_steady_state::Vector{Float64}
    endogenous_terminal_steady_state::Vector{Float64}
    endogenous_linear_trend::Vector{Float64}
    endogenous_quadratic_trend::Vector{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_terminal_steady_state::Vector{Float64}
    exogenous_linear_trend::Vector{Float64}
    exogenous_quadratic_trend::Vector{Float64}
    exogenous_det_steady_state::Vector{Float64}
    exogenous_det_terminal_steady_state::Vector{Float64}
    exogenous_det_linear_trend::Vector{Float64}
    exogenous_det_quadratic_trend::Vector{Float64}
    function Trends(ny::Int, nx::Int, nxd::Int)
        endogenous_steady_state = Vector{Float64}(undef, ny)
        endogenous_terminal_steady_state = Vector{Float64}(undef, 0)
        endogenous_linear_trend = zeros(ny)
        endogenous_quadratic_trend = Vector{Float64}(undef, ny)
        exogenous_steady_state = Vector{Float64}(undef, nx)
        exogenous_terminal_steady_state = Vector{Float64}(undef, 0)
        exogenous_linear_trend = Vector{Float64}(undef, nx)
        exogenous_quadratic_trend = Vector{Float64}(undef, nx)
        exogenous_det_steady_state = Vector{Float64}(undef, nxd)
        exogenous_det_terminal_steady_state = Vector{Float64}(undef, 0)
        exogenous_det_linear_trend = Vector{Float64}(undef, nxd)
        exogenous_det_quadratic_trend = Vector{Float64}(undef, nxd)
        new(
            endogenous_steady_state,
            endogenous_terminal_steady_state,
            endogenous_linear_trend,
            endogenous_quadratic_trend,
            exogenous_steady_state,
            exogenous_terminal_steady_state,
            exogenous_linear_trend,
            exogenous_quadratic_trend,
            exogenous_det_steady_state,
            exogenous_det_terminal_steady_state,
            exogenous_det_linear_trend,
            exogenous_det_quadratic_trend,
        )
    end
end

Base.show(io::IO, t::Trends) = show_field_value(t)

mutable struct EstimationResults
    posterior_mode
    posterior_mode_covariance
    posterior_mcmc_chains_nbr
    function EstimationResults()
        posterior_mode = []
        posterior_mode_covariance = Matrix(undef, 0, 0)
        posterior_mcmc_chains_nbr = 0
        new(posterior_mode, posterior_mode_covariance, posterior_mcmc_chains_nbr)
    end 
end 

mutable struct ModelResults
    irfs::Dict{Symbol, AxisArrayTable}
    trends::Trends
    stationary_variables::Vector{Bool}
    estimation::EstimationResults
    linearrationalexpectations::LinearRationalExpectationsResults
    simulations::Vector{Simulation}
    smoother
    solution_derivatives::Vector{Matrix{Float64}}
end

Base.show(io::IO, r::ModelResults) = show_field_value(r)

struct Results
    model_results::Vector{ModelResults}
end

struct EstimatedParameters
    index::Vector{Union{Int,Pair{Int,Int}}}
    initialvalue::Vector{Union{Float64, Missing}}
    ml_maximizer::Vector{Float64}
    name::Vector{Union{String,Pair{String,String}}}
    parametertype::Vector{EstimatedParameterType}
    posterior_mean::Vector{Float64}
    posterior_median::Vector{Float64}
    posterior_mode::Vector{Float64}
    posterior_sd::Vector{Float64}
    posterior_hpdi_lb::Vector{Float64}
    posterior_hpdi_ub::Vector{Float64}
    prior::Vector{Distribution}
    function EstimatedParameters()
        index = Vector{Union{Int,Pair{Int,Int}}}(undef, 0)
        initialvalue = Vector{Union{Float64, Missing}}(undef, 0)
        ml_maximizer = Vector{Float64}(undef, 0)
        name = Vector{Union{String,Pair{String,String}}}(undef, 0)
        parametertype = Vector{EstimatedParameterType}(undef, 0)
        posterior_mean = Vector{Float64}(undef, 0)
        posterior_median = Vector{Float64}(undef, 0)
        posterior_mode = Vector{Float64}(undef, 0)
        posterior_sd = Vector{Float64}(undef, 0)
        posterior_hpdi_lb = Vector{Float64}(undef, 0)
        posterior_hpdi_ub = Vector{Float64}(undef, 0)
        prior = Vector{Distribution}(undef, 0)
        new(
            index,
            initialvalue,
            ml_maximizer,
            name,
            parametertype,
            posterior_mean,
            posterior_median,
            posterior_mode,
            posterior_sd,
            posterior_hpdi_lb,
            posterior_hpdi_ub,
            prior,
        )
    end
end

import Base: length
Base.length(ep::EstimatedParameters) = length(ep.prior)

mutable struct Work
    analytical_steadystate_variables::Vector{Int}
    params::Vector{Float64}
    residuals::Vector{Float64}
    dynamic_variables::Vector{Float64}
    exogenous_variables::Vector{Float64}
    observed_variables::Vector{String}
    Sigma_m::Matrix{Float64}
    jacobian::Matrix{Float64}
    qr_jacobian::Matrix{Float64}
    model_has_trend::Vector{Bool}
    histval::Matrix{Union{Float64,Missing}}
    initval_endogenous::Matrix{Union{Float64,Missing}}
    initval_exogenous::Matrix{Union{Float64,Missing}}
    initval_exogenous_deterministic::Matrix{Union{Float64,Missing}}
    endval_endogenous::Matrix{Union{Float64,Missing}}
    endval_exogenous::Matrix{Union{Float64,Missing}}
    endval_exogenous_deterministic::Matrix{Union{Float64,Missing}}
    shocks::Vector{Float64}
    perfect_foresight_setup::Dict{String,Any}
    estimated_parameters::EstimatedParameters
    function Work(model, varobs)
        analytical_steadystate_variables = Int[]
        endo_nbr = model.endogenous_nbr
        exo_nbr = model.exogenous_nbr
        exo_det_nbr = model.exogenous_deterministic_nbr
        ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
        ncol1 = ncol + exo_nbr
        params = Vector{Float64}(undef, model.parameter_nbr)
        residuals = Vector{Float64}(undef, endo_nbr)
        dynamic_variables = Vector{Float64}(undef, ncol)
        # reserve enough space for a single period computation
        exogenous_variables = Vector{Float64}(undef, 3 * model.exogenous_nbr)
        observed_variables = varobs
        Sigma_m = zeros(length(varobs), length(varobs))
        #        jacobian = Matrix{Float64}(undef, model.endogenous_nbr, ncol1)
        #        qr_jacobian = Matrix{Float64}(undef, model.endogenous_nbr, ncol1)
        jacobian = Matrix{Float64}(undef, 0, 0)
        qr_jacobian = Matrix{Float64}(undef, 0, 0)
        model_has_trend = [false]
        histval = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        initval_endogenous = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        initval_exogenous = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        initval_exogenous_deterministic = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        endval_endogenous = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        endval_exogenous = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        endval_exogenous_deterministic = Matrix{Union{Float64,Missing}}(undef, 0, 0)
        # shocks
        shocks = Vector{Float64}(undef, 0)
        perfect_foresight_setup = Dict("periods" => 0, "datafile" => "")
        estimated_parameters = EstimatedParameters()
        new(
            analytical_steadystate_variables,
            params,
            residuals,
            dynamic_variables,
            exogenous_variables,
            observed_variables,
            Sigma_m,
            jacobian,
            qr_jacobian,
            model_has_trend,
            histval,
            initval_endogenous,
            initval_exogenous,
            initval_exogenous_deterministic,
            endval_endogenous,
            endval_exogenous,
            endval_exogenous_deterministic,
            shocks,
            perfect_foresight_setup,
            estimated_parameters,
        )
    end
end

Base.show(io::IO, w::Work) = show_field_value(w)

struct DynareSymbol
    longname::String
    texname::String
    symboltype::SymbolType
    orderintype::Integer
end

Base.getproperty(d::DynareSymbol, s::Symbol) = getfield(d, s)
Base.show(io::IO, ds::DynareSymbol) = show_field_value(ds)

const SymbolTable = Dict{String,DynareSymbol}

struct Context
    symboltable::SymbolTable
    models::Vector{Model}
    modfileinfo::ModFileInfo
    results::Results
    work::Work
    workspaces::Dict
end

# For now workspace assumes that only 1 model is present.
# In the future, if more models should exist, there should
# be a dict with the model id as key and each time a specific
# Ws for a specific model is needed that one is taken, or generated.
for ws_t in (:LinearCyclicReductionWs, :LinearGsSolverWs)
    @eval function workspace(::Type{LRE.$ws_t}, c::Context; model_id = 1)
        m = c.models[model_id]
        if !haskey(c.workspaces, LRE.$ws_t)
            c.workspaces[LRE.$ws_t] = LRE.$ws_t(m.ids)
        end
        return c.workspaces[LRE.$ws_t]
    end
end

# Probably should think of using an Algo type so that we can  dispatch rather than if elseif below
# IF more algos than 2 appear.
function workspace(::Type{LRE.LinearRationalExpectationsWs}, c::Context; model_id=1, algo="GS")
    m = c.models[model_id]
    if !haskey(c.workspaces, LRE.LinearRationalExpectationsWs)
        c.workspaces[LRE.LinearRationalExpectationsWs] = Dict{String, LRE.LinearRationalExpectationsWs}()
    end
    if !haskey(c.workspaces[LRE.LinearRationalExpectationsWs], algo)
        if algo == "GS"
            solver_ws = workspace(LRE.LinearGsSolverWs, c; model_id=model_id)
        elseif algo == "CR"
            solver_ws = workspace(LRE.LinearCyclicReductionWs, c; model_id=model_id)
        else
            error("Linear solving algorithm $algo not recognized, allowed: \"GS\", \"CR\"")
        end
        c.workspaces[LRE.LinearRationalExpectationsWs][algo] =
        LRE.LinearRationalExpectationsWs(solver_ws, m.ids)
    end
    return c.workspaces[LRE.LinearRationalExpectationsWs][algo]
end

function workspace(::Type{LRE.LyapdWs}, c::Context; model_id=1)
    m = c.models[model_id]
    if !haskey(c.workspaces, LRE.LyapdWs)
        c.workspaces[LRE.LyapdWs] = LRE.LyapdWs(m.n_bkwrd + m.n_both)
    end
    return c.workspaces[LRE.LyapdWs]
end
function workspace(::Type{LRE.VarianceWs}, c::Context; model_id=1, algo="GS")
    m = c.models[model_id]
    if !haskey(c.workspaces, LRE.VarianceWs)
        c.workspaces[LRE.VarianceWs] = Dict{String, LRE.VarianceWs}()
    end
    if !haskey(c.workspaces[LRE.VarianceWs], algo)
        lya = workspace(LRE.LyapdWs, c; model_id=model_id)
        lre = workspace(LRE.LinearRationalExpectationsWs, c; model_id=model_id, algo=algo)
        c.workspaces[LRE.VarianceWs][algo] = LRE.VarianceWs(m.endogenous_nbr, lya, m.exogenous_nbr, lre)
    end
    return c.workspaces[LRE.VarianceWs][algo]
end

"""
parameters obtained from modfile.json
"""
struct ModelInfo
    lead_lag_incidence::Vector{Any}
    nstatic::Int
    nfwrd::Int
    npred::Int
    nboth::Int
    nsfwrd::Int
    nspred::Int
    ndynamic::Int
    maximum_endo_lag::Int
    maximum_endo_lead::Int
    maximum_exo_lag::Int
    maximum_exo_lead::Int
    maximum_exo_det_lag::Int
    maximum_exo_det_lead::Int
    maximum_lag::Int
    maximum_lead::Int
    orig_maximum_endo_lag::Int
    orig_maximum_endo_lead::Int
    orig_maximum_exo_lag::Int
    orig_maximum_exo_lead::Int
    orig_maximum_exo_det_lag::Int
    orig_maximum_exo_det_lead::Int
    orig_maximum_lag::Int
    orig_maximum_lead::Int
    NNZDerivatives::Vector{Any}
end

Base.show(io::IO, m::ModelInfo) = show_field_value(m)

function show_field_value(s::Any)
    fn = fieldnames(typeof(s))
    for f in fn
        println("$f: $(getproperty(s, f))")
    end
end

#=
function Base.vcat(v1::Vector{T}, v2::Vector{T}) where {T}
    n1 = length(v1)
    n2 = length(v2)
    n = n1 + n2
    arr = Vector{T}(undef, n)
    unsafe_copyto!(arr, 1, v1, 1, n1)
    unsafe_copyto!(arr, n1 + 1, v2, 1, n2)
    return arr
end
=#

function Base.vcat(v1::Vector{T}, v2::Vector{T}, v3::Vector{T}) where {T}
    n1 = length(v1)
    n2 = length(v2)
    n3 = length(v3)
    n = n1 + n2 + n3
    arr = Vector{T}(undef, n)
    unsafe_copyto!(arr, 1, v1, 1, n1)
    unsafe_copyto!(arr, n1 + 1, v2, 1, n2)
    unsafe_copyto!(arr, n1 + n2 + 1, v3, 1, n3)
    return arr
end
#=
function load_dynare_function(modname::String, compileoption::Bool)#::Module
    @suppress begin
        fun = readlines(modname * ".jl")
        if fun[6] == "using StatsFuns"
            fun[6] = "using Dynare.StatsFuns"
        else
            insert!(fun, 6, "using Dynare.StatsFuns")
        end
        return eval(Meta.parse(join(fun, "\n")))
    end
end

function load_dynare_function2(modname::String)::Function
    fun = readlines(modname * ".jl")
    return (@RuntimeGeneratedFunction(Meta.parse(join(fun[3:(end-1)], "\n"))))
end
=#
