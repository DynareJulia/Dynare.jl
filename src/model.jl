using LinearAlgebra
using SparseArrays
using Suppressor

export get_de, get_abc, inverse_order_of_dynare_decision_rule

struct DynamicWs
    dynamic_variables::Vector{Float64}
    exogenous_variables::Vector{Float64}
    residuals::Vector{Float64}
    derivatives::Vector{SparseMatrixCSC}
    temporary_values::Vector{Float64}
    function DynamicWs(
        endogenous_nbr::Int64,
        exogenous_nbr::Int64,
        tmp_nbr::Int64,
        colptr::AbstractVector{Int64},
        rowval::AbstractVector{Int64}
    )
        dynamic_variables = Vector{Float64}(undef, 3*endogenous_nbr)
        exogenous_variables = Vector{Float64}(undef, exogenous_nbr)
        residuals = Vector{Float64}(undef, endogenous_nbr)
        derivatives = [SparseMatrixCSC(endogenous_nbr,
                                       3*endogenous_nbr + exogenous_nbr,
                                       colptr,
                                       rowval,
                                       similar(rowval, Float64))]
        temporary_values = Vector{Float64}(undef, tmp_nbr)
        new(
            dynamic_variables,
            exogenous_variables,
            residuals,
            derivatives,
            temporary_values,
        )
    end
end

function DynamicWs(context::Context)
    m = context.models[1]
    tmp_nbr = sum(m.dynamic_tmp_nbr[1:2])
    return DynamicWs(m.endogenous_nbr,
                     m.exogenous_nbr,
                     tmp_nbr,
                     m.dynamic_g1_sparse_colptr,
                     m.dynamic_g1_sparse_rowval)
end

struct StaticWs
    residuals::Vector{Float64}
    derivatives::Vector{SparseMatrixCSC}
    temporary_values::Vector{Float64}
    function StaticWs(
        endogenous_nbr::Int64,
        tmp_nbr::Int64,
        colptr::AbstractVector{Int64},
        rowval::AbstractVector{Int64}
    )
        residuals = Vector{Float64}(undef, endogenous_nbr)
        derivatives = [SparseMatrixCSC(endogenous_nbr,
                                       endogenous_nbr,
                                       colptr,
                                       rowval,
                                       similar(rowval, Float64))]
        temporary_values = Vector{Float64}(undef, tmp_nbr)
        new(residuals, derivatives, temporary_values)
    end
end

function StaticWs(context::Context)
    m = context.models[1]
    tmp_nbr = sum(m.static_tmp_nbr[1:2])
    return StaticWs(m.endogenous_nbr,
                    tmp_nbr,
                    m.static_g1_sparse_colptr,
                    m.static_g1_sparse_rowval)
end

function get_exogenous_matrix(x::Vector{Float64}, exogenous_nbr::Int64)
    @debug "any(isnan.(x))=$(any(isnan.(x))) "
    x1 = reshape(x, Int(length(x) / exogenous_nbr), exogenous_nbr)
    @debug "any(isnan.(x1))=$(any(isnan.(x1))) "
    return x1
end

function get_static_residuals!(
    ws::StaticWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
#    df::DynareFunctions,
)
    DFunctions.static!(
        ws.temporary_values,
        ws.residuals,
        endogenous,
        exogenous,
        params,
    )
    return ws.residuals
end

"""
`get_dynamic_jacobian!`(ws::DynamicWs, params::Vector{Float64}, endogenous::AbstractVector{Float64}, exogenous::Vector{Float64}, m::Model, period::Int64)

sets the dynamic Jacobian matrix ``work.jacobian``, evaluated at ``endogenous`` and ``exogenous`` values, identical for all leads and lags
"""
function get_dynamic_jacobian!(
    ws::DynamicWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::Vector{Float64},
    m::Model,
#    df::DynareFunctions,
    period::Int64,
)
    DFunctions.dynamic!(
        ws.temporary_values,
        ws.residuals,
        ws.derivatives[1],
        endogenous,
        exogenous,
        params,
        steadystate,
    )
    return ws.derivatives[1]
end

"""
`get_static_jacobian!`(ws::StaticWs, params::Vector{Float64}, endogenous::AbstractVector{Float64}, exogenous::Vector{Float64})

sets the static Jacobian matrix ``work.jacobian``, evaluated at ``endogenous`` and ``exogenous`` values
"""
function get_static_jacobian!(
    ws::StaticWs,
    params::Vector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    m::Model,
#    df::DynareFunctions,
)
    @debug "any(isnan.(exognous))=$(any(isnan.(exogenous)))"
    DFunctions.static!(
        ws.temporary_values,
        ws.residuals,
        ws.derivatives[1],
        endogenous,
        exogenous,
        params,
    )
    return ws.derivatives[1]
end

function get_abc!(
    a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64},
    jacobian::AbstractMatrix{Float64},
    m::Model,
)
    i_rows = (m.n_static+1):m.endogenous_nbr
    fill!(a, 0.0)
    fill!(b, 0.0)
    fill!(c, 0.0)
    ws.a[:, ws.forward_indices_d] .=
        view(jacobian, i_rows, ws.backward_nbr .+ ws.current_nbr .+ (1:ws.forward_nbr))
    ws.b[:, ws.current_dynamic_indices_d] .=
        view(jacobian, i_rows, ws.backward_nbr .+ ws.current_dynamic_indices)
    ws.c[:, ws.backward_indices_d] .= view(jacobian, i_rows, 1:ws.backward_nbr)
end

function get_de!(ws::LinearRationalExpectationsWs, jacobian::AbstractMatrix{Float64})
    n1 = ws.backward_nbr + ws.forward_nbr - ws.both_nbr
    fill!(ws.d, 0.0)
    fill!(ws.e, 0.0)
    i_rows = (ws.static_nbr+1):ws.endogenous_nbr
    ws.d[1:n1, ws.icolsD] .= jacobian[i_rows, ws.jcolsD]
    ws.e[1:n1, ws.icolsE] .= -jacobian[i_rows, ws.jcolsE]
    u = Matrix{Float64}(I, ws.both_nbr, ws.both_nbr)
    i_rows = n1 .+ (1:ws.both_nbr)
    ws.d[i_rows, ws.colsUD] .= u
    ws.e[i_rows, ws.colsUE] .= u
end



