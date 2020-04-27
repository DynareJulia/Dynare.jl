using LinearAlgebra
export Model, get_de, get_abc, inverse_order_of_dynare_decision_rule

struct Model
    endogenous_nbr
    exogenous_nbr
    lagged_exogenous_nbr
    exogenous_deterministic_nbr
    parameter_nbr
    lead_lag_incidence
    n_static
    n_fwrd
    n_bkwrd
    n_both
    n_states
    DErows1
    DErows2
    n_dyn
    i_static
    i_dyn::Array{Int64,1}
    i_bkwrd
    i_bkwrd_b
    i_bkwrd_ns
    i_fwrd
    i_fwrd_b
    i_fwrd_ns
    i_both
    p_static
    p_bkwrd
    p_bkwrd_b
    p_fwrd
    p_fwrd_b
    p_both_b
    p_both_f
    i_current
    p_current
    n_current
    i_current_ns
    p_current_ns
    n_current_ns
    icolsD
    jcolsD
    icolsE
    jcolsE
    colsUD
    colsUE
    i_cur_fwrd
    n_cur_fwrd
    p_cur_fwrd
    i_cur_bkwrd
    n_cur_bkwrd
    p_cur_bkwrd
    i_cur_both
    n_cur_both
    p_cur_both
    gx_rows
    hx_rows
    i_current_exogenous
    i_lagged_exogenous
    serially_correlated_exogenous
    dynamic!
    static!
    steady_state!
end

function Model(modfilename, endo_nbr, lead_lag_incidence, exogenous_nbr, lagged_exogenous_nbr, exogenous_deterministic_nbr, parameter_nbr)
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
    icolsD = [1:n_cur_bkwrd; n_bkwrd + n_both .+ (1:(n_fwrd+n_both))]
    jcolsD = [p_cur_bkwrd; p_fwrd; p_both_f]
    # derivatives of current values of variables that are both
    # forward and backward are included in the E matrix
    icolsE = [1:(n_bkwrd + n_both); n_bkwrd + n_both .+ (1:(n_fwrd+n_both))]
    jcolsE = [p_bkwrd; p_both_b; p_cur_fwrd; p_cur_both]
    colsUD = n_bkwrd .+ (1:n_both)
    colsUE = n_both + n_fwrd .+ colsUD
    n_dyn = endo_nbr - n_static + n_both
    DErows1 = 1:(n_dyn-n_both)
    DErows2 = (n_dyn-n_both) .+ (1:n_both)
    gx_rows = n_bkwrd .+ (1:(n_fwrd+n_both))
    hx_rows = 1:(n_bkwrd + n_both)
    i_current_exogenous = maximum(lead_lag_incidence) .+ (1:exogenous_nbr)
    i_lagged_exogenous = 0:-1
    serially_correlated_exogenous = false
    dynamic! = load_dynare_function(modfilename*"Dynamic.jl")
    static! = load_dynare_function(modfilename*"Static.jl")
    if isfile(modfilename*"SteadyState2.jl")
        steady_state! = load_dynare_function(modfilename*"SteadyState2.jl")
    else
        steady_state! = nothing
    end
    Model(endo_nbr, exogenous_nbr, lagged_exogenous_nbr, exogenous_deterministic_nbr,
          parameter_nbr, lead_lag_incidence, n_static, n_fwrd, n_bkwrd, n_both,
          n_states, DErows1, DErows2, n_dyn, i_static, i_dyn, i_bkwrd,
          i_bkwrd_b, i_bkwrd_ns, i_fwrd, i_fwrd_b, i_fwrd_ns, i_both,
          p_static, p_bkwrd, p_bkwrd_b, p_fwrd, p_fwrd_b, p_both_b,
          p_both_f, i_current, p_current, n_current, i_current_ns,
          p_current_ns, n_current_ns, icolsD, jcolsD, icolsE, jcolsE,
          colsUD, colsUE, i_cur_fwrd, n_cur_fwrd, p_cur_fwrd,
          i_cur_bkwrd, n_cur_bkwrd, p_cur_bkwrd, i_cur_both,
          n_cur_both, p_cur_both, gx_rows, hx_rows,
          i_current_exogenous, i_lagged_exogenous,
          serially_correlated_exogenous, dynamic!, static!,
          steady_state!)   
end

function inverse_order_of_dynare_decision_rule(m::Model)
    inverse_order_var = Vector{Int64}(undef, m.endo_nbr)
    for i = 1:m.n_static
        inverse_order_var[m.i_static[i]] = i
    end

    offset = m.n_static
    for i = 1:m.n_bkwrd
        inverse_order_var[m.i_bkwrd[i]] = i + offset
    end

    offset += m.n_bkwrd
    for i = 1:m.n_both
        inverse_order_var[m.i_both[i]] = i + offset
    end

    offset += m.n_both
    for i = 1:m.n_fwrd
        inverse_order_var[m.i_fwrd[i]] = i + offset
    end

    inverse_order_states = sortperm(cat(m.i_bkwrd,m.i_both;dims=1))

    (inverse_order_var, inverse_order_states)
end

function load_dynare_function(filename)
    file = readlines(filename)
    # drop using Utils
    file[6] = "using Dynare: get_power_deriv"
    return eval(Meta.parse(join(file, "\n")))
end

