import Base
using LinearRationalExpectations
using TimeDataFrames

export Context, DynareSymbol, Model, ModelResults, Results, Simulation, SymbolType, Work, Trends

struct Model
    endogenous_nbr::Int64
    exogenous_nbr::Int64
    lagged_exogenous_nbr::Int64
    exogenous_deterministic_nbr::Int64
    parameter_nbr::Int64
    lead_lag_incidence::Matrix{Int64}
    n_static::Int64
    n_fwrd::Int64
    n_bkwrd::Int64
    n_both::Int64
    n_states::Int64
    DErows1::Vector{Int64}
    DErows2::Vector{Int64}
    n_dyn::Int64
    i_static::Vector{Int64}
    i_dyn::Array{Int64,1}
    i_bkwrd::Vector{Int64}
    i_bkwrd_b::Vector{Int64}
    i_bkwrd_ns::Vector{Int64}
    i_fwrd::Vector{Int64}
    i_fwrd_b::Vector{Int64}
    i_fwrd_ns::Vector{Int64}
    i_both::Vector{Int64}
    i_non_states::Vector{Int64}
    p_static::Vector{Int64}
    p_bkwrd::Vector{Int64}
    p_bkwrd_b::Vector{Int64}
    p_fwrd::Vector{Int64}
    p_fwrd_b::Vector{Int64}
    p_both_b::Vector{Int64}
    p_both_f::Vector{Int64}
    i_current::Vector{Int64}
    p_current::Vector{Int64}
    n_current::Int64
    i_current_ns::Vector{Int64}
    p_current_ns::Vector{Int64}
    n_current_ns::Int64
    icolsD::Vector{Int64}
    jcolsD::Vector{Int64}
    icolsE::Vector{Int64}
    jcolsE::Vector{Int64}
    colsUD::Vector{Int64}
    colsUE::Vector{Int64}
    i_cur_fwrd::Vector{Int64}
    n_cur_fwrd::Int64
    p_cur_fwrd::Vector{Int64}
    i_cur_bkwrd::Vector{Int64}
    n_cur_bkwrd::Int64
    p_cur_bkwrd::Vector{Int64}
    i_cur_both::Vector{Int64}
    n_cur_both::Int64
    p_cur_both::Vector{Int64}
    gx_rows::Vector{Int64}
    hx_rows::Vector{Int64}
    i_current_exogenous::Vector{Int64}
    i_lagged_exogenous::Vector{Int64}
    serially_correlated_exogenous::Vector{Int64}
    Sigma_e::Matrix{Float64}
    maximum_endo_lag::Int64
    maximum_endo_lead::Int64
    maximum_exo_lag::Int64
    maximum_exo_lead::Int64
    maximum_exo_det_lag::Int64
    maximum_exo_det_lead::Int64
    maximum_lag::Int64
    maximum_lead::Int64
    orig_maximum_endo_lag::Int64
    orig_maximum_endo_lead::Int64
    orig_maximum_exo_lag::Int64
    orig_maximum_exo_lead::Int64
    orig_maximum_exo_det_lag::Int64
    orig_maximum_exo_det_lead::Int64
    orig_maximum_lag::Int64
    orig_maximum_lead::Int64
    dynamic_indices::Vector{Int64}
    current_dynamic_indices::Vector{Int64}
    forward_indices_d::Vector{Int64}
    backward_indices_d::Vector{Int64}
    current_dynamic_indices_d::Vector{Int64}
    exogenous_indices::Vector{Int64}
    NNZDerivatives::Vector{Int64}
    dynamic!::Module
    static!::Module
    set_auxiliary_variables!::Module
    set_dynamic_auxiliary_variables!::Module
    steady_state!::Module
end

function Model(modfilename::String, endogenous_nbr::Int64,
               lead_lag_incidence::Matrix{Int64},
               exogenous_nbr::Int64, lagged_exogenous_nbr::Int64, exogenous_deterministic_nbr::Int64,
               parameter_nbr::Int64, maximum_endo_lag::Int64, maximum_endo_lead::Int64,
               maximum_exo_lag::Int64, maximum_exo_lead::Int64, maximum_exo_det_lag::Int64,
               maximum_exo_det_lead::Int64, maximum_lag, maximum_lead::Int64,
               orig_maximum_endo_lag::Int64, orig_maximum_endo_lead::Int64,
               orig_maximum_exo_lag::Int64, orig_maximum_exo_lead::Int64,
               orig_maximum_exo_det_lag::Int64, orig_maximum_exo_det_lead::Int64,
               orig_maximum_lag::Int64, orig_maximum_lead::Int64, NNZDerivatives::Vector{Int64}, compileoption::Bool)
    i_static = findall((lead_lag_incidence[1,:] .== 0) .& (lead_lag_incidence[3,:] .== 0))
    p_static = lead_lag_incidence[2,i_static]
    i_dyn = findall((lead_lag_incidence[1,:] .> 0) .| (lead_lag_incidence[3,:] .> 0))
    n_static = length(i_static)
    i_bkwrd = findall((lead_lag_incidence[1,:] .> 0) .& (lead_lag_incidence[3,:] .== 0))
    i_bkwrd_b = findall((lead_lag_incidence[1,:] .> 0))
    i_bkwrd_ns = findall(lead_lag_incidence[1,i_dyn] .> 0)
    p_bkwrd = lead_lag_incidence[1,i_bkwrd]
    p_bkwrd_b = lead_lag_incidence[1,i_bkwrd_b]
    n_bkwrd = length(i_bkwrd)
    i_fwrd = findall((lead_lag_incidence[3,:] .> 0) .& (lead_lag_incidence[1,:] .== 0)) 
    i_fwrd_b = findall((lead_lag_incidence[3,:] .> 0)) 
    i_fwrd_ns = findall(lead_lag_incidence[3,i_dyn] .> 0)
    p_fwrd = lead_lag_incidence[3,i_fwrd]
    p_fwrd_b = lead_lag_incidence[3,i_dyn[i_fwrd_ns]]
    n_fwrd = length(i_fwrd)
    i_both = findall((lead_lag_incidence[1,:] .> 0) .& (lead_lag_incidence[3,:] .> 0))
    i_non_states = union(i_fwrd, i_static)
    p_both_b = lead_lag_incidence[1,i_both]
    p_both_f = lead_lag_incidence[3,i_both]
    n_both = length(i_both)
    n_states = n_bkwrd + n_both
    i_current = findall(lead_lag_incidence[2,:] .> 0 )
    p_current = lead_lag_incidence[2,i_current]
    n_current = count(i->(i > 0),lead_lag_incidence[2,:])
    i_current_ns = findall(lead_lag_incidence[2,i_dyn] .> 0 )
    p_current_ns = lead_lag_incidence[2,i_dyn[i_current_ns]]
    n_current_ns = count(i->(i > 0),lead_lag_incidence[2,i_dyn])
    i_cur_fwrd = findall(lead_lag_incidence[2,i_fwrd] .> 0)
    n_cur_fwrd = length(i_cur_fwrd)
    p_cur_fwrd = lead_lag_incidence[2,i_fwrd[i_cur_fwrd]]
    i_cur_bkwrd = findall(lead_lag_incidence[2,i_bkwrd] .> 0)
    n_cur_bkwrd = length(i_cur_bkwrd)
    p_cur_bkwrd = lead_lag_incidence[2,i_bkwrd[i_cur_bkwrd]]
    i_cur_both = findall(lead_lag_incidence[2,i_both] .> 0)
    n_cur_both = length(i_cur_both)
    p_cur_both = lead_lag_incidence[2,i_both[i_cur_both]]
    i_backward_in_current = findall(in(i_current),i_bkwrd)
    icolsD = [i_backward_in_current;
              n_bkwrd + n_both .+ collect(1:(n_fwrd+n_both))]
    jcolsD = [p_cur_bkwrd; p_fwrd; p_both_f]
    # derivatives of current values of variables that are both
    # forward and backward are included in the E matrix
    icolsE = [collect(1:(n_bkwrd + n_both)); n_bkwrd + n_both .+ collect(1:(n_fwrd+n_both))]
    jcolsE = [p_bkwrd; p_both_b; p_cur_fwrd; p_cur_both]
    colsUD = n_bkwrd .+ collect(1:n_both)
    colsUE = n_both + n_fwrd .+ colsUD
    n_dyn = endogenous_nbr - n_static + n_both
    DErows1 = collect(1:(n_dyn-n_both))
    DErows2 = (n_dyn-n_both) .+ collect(1:n_both)
    gx_rows = n_bkwrd .+ collect(1:(n_fwrd+n_both))
    hx_rows = collect(1:(n_bkwrd + n_both))
    i_current_exogenous = maximum(lead_lag_incidence) .+ (1:exogenous_nbr)
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
    current_number = n_current
    dynamic_indices = setdiff(collect(1:endogenous_nbr), static_indices)
    current_dynamic_indices = setdiff(current_indices, static_indices)
    purely_forward_indices = setdiff(forward_indices, both_indices)
    forward_indices_d = findall(in(forward_indices), dynamic_indices)
    backward_indices_d = findall(in(backward_indices), dynamic_indices)
    current_dynamic_indices_d = findall(in(current_dynamic_indices), dynamic_indices)
    exogenous_indices = (backward_number + current_number
                         + forward_number .+ collect(1:exogenous_nbr))
    
    dynamic! = load_dynare_function(modfilename*"Dynamic",
                                    compileoption)
    static! = load_dynare_function(modfilename*"Static",
                                   compileoption)
    if isfile(modfilename*"DynamicSetAuxiliarySeries.jl")
        set_dynamic_auxiliary_variables! =
            load_dynare_function(modfilename*"DynamicSetAuxiliarySeries",
                                 compileoption)
    else
        set_dynamic_auxiliary_variables! = Module()
    end
    if isfile(modfilename*"SetAuxiliaryVariables.jl")
        set_auxiliary_variables! =
            load_dynare_function(modfilename*"SetAuxiliaryVariables",
                                 compileoption)
    else
        set_auxiliary_variables! = Module()
    end
    if isfile(modfilename*"SteadyState2.jl")
        steady_state! = load_dynare_function(modfilename*"SteadyState2",
                                             compileoption)
    else
        steady_state! = Module()
    end
    Model(endogenous_nbr, exogenous_nbr, lagged_exogenous_nbr,
          exogenous_deterministic_nbr, parameter_nbr,
          lead_lag_incidence, n_static, n_fwrd, n_bkwrd, n_both,
          n_states, DErows1, DErows2, n_dyn, i_static, i_dyn, i_bkwrd,
          i_bkwrd_b, i_bkwrd_ns, i_fwrd, i_fwrd_b, i_fwrd_ns, i_both,
          i_non_states, p_static, p_bkwrd, p_bkwrd_b, p_fwrd,
          p_fwrd_b, p_both_b, p_both_f, i_current, p_current,
          n_current, i_current_ns, p_current_ns, n_current_ns, icolsD,
          jcolsD, icolsE, jcolsE, colsUD, colsUE, i_cur_fwrd,
          n_cur_fwrd, p_cur_fwrd, i_cur_bkwrd, n_cur_bkwrd,
          p_cur_bkwrd, i_cur_both, n_cur_both, p_cur_both, gx_rows,
          hx_rows, i_current_exogenous, i_lagged_exogenous,
          serially_correlated_exogenous, Sigma_e, maximum_endo_lag,
          maximum_endo_lead, maximum_exo_lag, maximum_exo_lead,
          maximum_exo_det_lag, maximum_exo_det_lead, maximum_lag,
          maximum_lead, orig_maximum_endo_lag, orig_maximum_endo_lead,
          orig_maximum_exo_lag, orig_maximum_exo_lead,
          orig_maximum_exo_det_lag, orig_maximum_exo_det_lead,
          orig_maximum_lag, orig_maximum_lead, dynamic_indices,
          current_dynamic_indices, forward_indices_d,
          backward_indices_d, current_dynamic_indices_d,
          exogenous_indices, NNZDerivatives, dynamic!, static!,
          set_auxiliary_variables!, set_dynamic_auxiliary_variables!,
          steady_state!)
end

Base.show(io::IO, m::Model) = show_field_value(m)

struct Simulation
    name::String
    statement::String
    options::Dict{String, Any}
    data::TimeDataFrame
end

Base.show(io::IO, s::Simulation) = show_field_value(s)

struct Trends
    endogenous_steady_state::Vector{Float64}
    endogenous_linear_trend::Vector{Float64}
    endogenous_quadratic_trend::Vector{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_linexbar_trend::Vector{Float64}
    exogenous_quadratic_trend::Vector{Float64}
    exogenous_det_steady_state::Vector{Float64}
    exogenous_det_linear_trend::Vector{Float64}
    exogenous_det_quadratic_trend::Vector{Float64}
    function Trends(ny::Int64, nx::Int64, nxd::Int64)
        endogenous_steady_state = Vector{Float64}(undef, ny)
        endogenous_linear_trend = Vector{Float64}(undef, ny)
        endogenous_quadratic_trend = Vector{Float64}(undef, ny)
        exogenous_steady_state = Vector{Float64}(undef, nx)
        exogenous_linear_trend = Vector{Float64}(undef, nx)
        exogenous_quadratic_trend = Vector{Float64}(undef, nx)
        exogenous_det_steady_state = Vector{Float64}(undef, nxd)
        exogenous_det_linear_trend = Vector{Float64}(undef, nxd)
        exogenous_det_quadratic_trend = Vector{Float64}(undef, nxd)
        new(endogenous_steady_state, endogenous_linear_trend,
            endogenous_quadratic_trend, exogenous_steady_state,
            exogenous_linear_trend, exogenous_quadratic_trend,
            exogenous_det_steady_state, exogenous_det_linear_trend,
            exogenous_det_quadratic_trend)
    end
end

Base.show(io::IO, t::Trends) = show_field_value(t)

mutable struct ModelResults
    endogenous_steady_state::Vector{Float64}
    trends::Trends
    endogenous_variance::Matrix{Float64}
    stationary_variables::Vector{Bool}
    exogenous_steady_state::Vector{Float64}
    exogenous_deterministic_steady_state::Vector{Float64}
    linearrationalexpectations::LinearRationalExpectationsResults
    simulations::Vector{Simulation}
    smoother::Dict{String, Any}
end

Base.show(io::IO, r::ModelResults) = show_field_value(r)

struct Results
    model_results::Vector{ModelResults}
end

mutable struct Work
    params::Vector{Float64}
    residuals::Vector{Float64}
    temporary_values::Vector{Float64}
    dynamic_variables::Vector{Float64}
    exogenous_variables::Matrix{Float64}
    observed_variables::Vector{String}
    jacobian::Matrix{Float64}
    qr_jacobian::Matrix{Float64}
    model_has_trend::Bool
    histval::Matrix{Float64}
end

Base.show(io::IO, w::Work) = show_field_value(w)

@enum SymbolType Endogenous Exogenous ExogenousDeterministic Parameter DynareFunction

struct DynareSymbol
    longname::String
    texname::String
    symboltype::SymbolType
    orderintype::Integer
end

Base.show(io::IO, ds::DynareSymbol) = show_field_value(ds)    

const SymbolTable = Dict{String, DynareSymbol}

mutable struct Context
    symboltable::SymbolTable
    models::Vector{Model}
    options::Dict
    results::Results
    work::Work
end

struct ModelInfo
    lead_lag_incidence::Matrix{Int64}
    nstatic::Int64
    nfwrd::Int64
    npred::Int64
    nboth::Int64
    nsfwrd::Int64
    nspred::Int64
    ndynamic::Int64
    maximum_endo_lag::Int64
    maximum_endo_lead::Int64
    maximum_exo_lag::Int64
    maximum_exo_lead::Int64
    maximum_exo_det_lag::Int64
    maximum_exo_det_lead::Int64
    maximum_lag::Int64
    maximum_lead::Int64
    orig_maximum_endo_lag::Int64
    orig_maximum_endo_lead::Int64
    orig_maximum_exo_lag::Int64
    orig_maximum_exo_lead::Int64
    orig_maximum_exo_det_lag::Int64
    orig_maximum_exo_det_lead::Int64
    orig_maximum_lag::Int64
    orig_maximum_lead::Int64
    NNZDerivatives::Vector{Int64}
end

Base.show(io::IO, m::ModelInfo) = show_field_value(m)

function show_field_value(s::Any)
    fn = fieldnames(typeof(s))
    for f in fn
        println("$f: $(getproperty(s, f))")
    end
end
