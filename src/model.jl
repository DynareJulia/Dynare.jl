using LinearAlgebra
export get_de, get_abc, inverse_order_of_dynare_decision_rule

function inverse_order_of_dynare_decision_rule(m::Model)
    inverse_order_var = Vector{Int64}(undef, m.endogenous_nbr)
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

function load_dynare_function(filename::String)
    file = readlines(filename)
    # drop using Utils
    file[6] = "using Dynare: get_power_deriv"
    return eval(Meta.parse(join(file, "\n")))
end

"""
get_dynamic_endogenous_variables!(y::Vector{Float64}, data::Vector{Float64}, lli::Matrix{Int64})

sets the vector of dynamic variables ``y``, evaluated at the same values 
for all leads and lags and taken in ``data`` vector 
"""
function get_dynamic_endogenous_variables!(y::Vector{Float64}, data::Vector{Float64}, lli::Matrix{Int64})
    for i = 1:size(lli,2)
        value = data[i]
        for j = 1:size(lli,1)
            k = lli[j, i]
            if k > 0
                y[k] = value
            end
        end
    end
end

"""
get_dynamic_endogenous_variables!(y::Vector{Float64}, data::Matrix{Float64}, lli::Matrix{Int64})

sets the vector of dynamic variables ``y`` with values in as many rows of ``data`` matrix
as there are leads and lags in the model. ``period`` is the current period.
"""
function get_dynamic_endogenous_variables!(y::Vector{Float64}, data::Matrix{Float64}, lli::Matrix{Int64}, m::Model, period::Int64)
    for i = 1:size(lli,2)
        p = period - m.maximum_lag - 1
        for j = 1:size(lli,1)
            k = lli[j, i]
            if k > 0
                y[k] = data[p + j, i]
            end
        end
    end
end

"""
get_jacobian!(work::Work, endogenous::Vector{Float64}, exogenous::Vector{Float64}, m::Model, period::Int64)

returns sets the Jacobian matrix ``work.jacobian``, evaluated at ``endogenous`` and ``exogenous`` values, identical for all leads and lags
"""
function get_jacobian!(work::Work, endogenous::Vector{Float64}, exogenous::Vector{Float64},
                       steadystate::Vector{Float64}, m::Model, period::Int64)
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables!(work.dynamic_variables, endogenous, lli)
    work.exogenous_variables = repeat(transpose(exogenous), size(lli, 1), 1)
    Base.invokelatest(m.dynamic!.dynamic!,
                      work.temporary_values,
                      work.residuals,
                      work.jacobian,
                      work.dynamic_variables,
                      work.exogenous_variables,
                      work.params,
                      steadystate,
                      period)  
end

"""
get_jacobian!(work::Work, endogenous::Matrix{Float64}, exogenous::Matrix{Float64}, m::Model, period::Int64)

returns sets the Jacobian matrix ``work.jacobian``, evaluated with ``endogenous`` and ``exogenous`` values taken
around ``period`` 
"""
function get_jacobian!(work::Work, endogenous::Matrix{Float64}, exogenous::Matrix{Float64}, steadystate::Vector{Float64}, m::Model, period::Int64)
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables!(work.dynamic_variables, endogenous, lli, m, period)
    Base.invokelatest(m.dynamic!.dynamic!,
                      work.temporary_values,
                      work.residuals,
                      work.jacobian,
                      work.dynamic_variables,
                      exogenous,
                      work.params,
                      steadystate,
                      period)  
end

function get_abc!(a::AbstractMatrix{Float64},
                  b::AbstractMatrix{Float64},
                  c::AbstractMatrix{Float64},
                  jacobian::AbstractMatrix{Float64},
                  m::Model)
    i_rows = (m.n_static + 1):m.endogenous_nbr
    fill!(a, 0.0)
    fill!(b, 0.0)
    fill!(c, 0.0)
    ws.a[:, ws.forward_indices_d] .= view(jacobian, i_rows, ws.backward_nbr .+ ws.current_nbr .+ (1:ws.forward_nbr))
    ws.b[:, ws.current_dynamic_indices_d] .= view(jacobian, i_rows, ws.backward_nbr .+ ws.current_dynamic_indices)
    ws.c[:, ws.backward_indices_d] .= view(jacobian, i_rows, 1:ws.backward_nbr)
end

function get_de!(ws, jacobian::AbstractMatrix{Float64})
    n1 = ws.backward_nbr + ws.forward_nbr - ws.both_nbr
    fill!(ws.d, 0.0)
    fill!(ws.e, 0.0)
    i_rows = (ws.static_nbr + 1):ws.endogenous_nbr
    ws.d[1:n1, ws.icolsD] .= jacobian[i_rows, ws.jcolsD]
    ws.e[1:n1, ws.icolsE] .= -jacobian[i_rows, ws.jcolsE]
    u = Matrix{Float64}(I, ws.both_nbr, ws.both_nbr)                                    
    i_rows = n1 .+ (1:ws.both_nbr)
    ws.d[i_rows, ws.colsUD] .= u
    ws.e[i_rows, ws.colsUE] .= u
end



