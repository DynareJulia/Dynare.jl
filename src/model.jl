using LinearAlgebra
export Model, get_de, get_abc, inverse_order_of_dynare_decision_rule

struct Model
    endo_nbr
    current_exogenous_nbr
    lagged_exogenous_nbr
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
end

function Model(endo_nbr, lead_lag_incidence, current_exogenous_nbr, lagged_exogenous_nbr)
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
    icolsD = [1:n_cur_bkwrd; n_bkwrd+n_both .+ (1:(n_fwrd+n_both))]
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
    i_current_exogenous = maximum(lead_lag_incidence) .+ (1:current_exogenous_nbr)
    i_lagged_exogenous = 0:-1
    serially_correlated_exogenous = false
    Model(endo_nbr, current_exogenous_nbr, lagged_exogenous_nbr,
          lead_lag_incidence, n_static, n_fwrd, n_bkwrd, n_both,
          n_states, DErows1, DErows2, n_dyn, i_static, i_dyn, i_bkwrd,
          i_bkwrd_b, i_bkwrd_ns, i_fwrd, i_fwrd_b, i_fwrd_ns, i_both,
          p_static, p_bkwrd, p_bkwrd_b, p_fwrd, p_fwrd_b, p_both_b,
          p_both_f, i_current, p_current, n_current, i_current_ns,
          p_current_ns, n_current_ns, icolsD, jcolsD, icolsE, jcolsE,
          colsUD, colsUE, i_cur_fwrd, n_cur_fwrd, p_cur_fwrd,
          i_cur_bkwrd, n_cur_bkwrd, p_cur_bkwrd, i_cur_both,
          n_cur_both, p_cur_both, gx_rows, hx_rows,
          i_current_exogenous, i_lagged_exogenous,
          serially_correlated_exogenous)   
end

Model(endo_nbr, lli, current_exogenous_nbr) = Model(endo_nbr, lli, current_exogenous_nbr, 0)
    
function get_de(jacobian,model)
    n1 = size(model.DErows1,1)
    n2 = model.n_dyn - n1;
    d = zeros(model.n_dyn,model.n_dyn)
    e = zeros(model.n_dyn,model.n_dyn)
    d[1:n1,model.icolsD] = jacobian[:,model.jcolsD]
    e[1:n1,model.icolsE] = -jacobian[:,model.jcolsE]
    u = Matrix{Float64}(I, n2, n2)                                    
    d[model.DErows2,model.colsUD] = u
    e[model.DErows2,model.colsUE] = u
    return d, e
end

function get_abc(model::Model,jacobian::Array{Float64})
    i_rows = model.n_static+1:model.endo_nbr
    n = length(i_rows)
    a = zeros(Float64,n,n)
    b = zeros(Float64,n,n)
    c = zeros(Float64,n,n)
    a[:,model.i_bkwrd_ns] = view(jacobian,i_rows,model.p_bkwrd_b)
    b[:,model.i_current_ns]  = view(jacobian,i_rows,model.p_current_ns)
    c[:,model.i_fwrd_ns]  = view(jacobian,i_rows,model.p_fwrd_b)
    return a, b, c
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


