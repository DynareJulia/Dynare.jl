module DFunctions

using RuntimeGeneratedFunctions
using SparseArrays
using StatsFuns
using AxisArrayTables

RuntimeGeneratedFunctions.init(@__MODULE__)

function load_model_functions(modelname::String)
    function_root = "$(modelname)/model/julia/"
    global SparseDynamicG1! = load_dynare_function("$(function_root)SparseDynamicG1!.jl")
    global SparseDynamicG1TT! =
        load_dynare_function("$(function_root)SparseDynamicG1TT!.jl")
    global SparseDynamicG2! = load_dynare_function("$(function_root)SparseDynamicG2!.jl")
    global SparseDynamicG2TT! =
        load_dynare_function("$(function_root)SparseDynamicG2TT!.jl")
    global SparseDynamicG3! = load_dynare_function("$(function_root)SparseDynamicG3!.jl")
    global SparseDynamicG3TT! =
        load_dynare_function("$(function_root)SparseDynamicG3TT!.jl")
    global SparseDynamicResid! =
        load_dynare_function("$(function_root)SparseDynamicResid!.jl")
    global SparseDynamicResidTT! =
        load_dynare_function("$(function_root)SparseDynamicResidTT!.jl")
    global SparseStaticG1! = load_dynare_function("$(function_root)SparseStaticG1!.jl")
    global SparseStaticG1TT! = load_dynare_function("$(function_root)SparseStaticG1TT!.jl")
    global SparseStaticG2! = load_dynare_function("$(function_root)SparseStaticG2!.jl")
    global SparseStaticG2TT! = load_dynare_function("$(function_root)SparseStaticG2TT!.jl")
    global SparseStaticG3! = load_dynare_function("$(function_root)SparseStaticG3!.jl")
    global SparseStaticG3TT! = load_dynare_function("$(function_root)SparseStaticG3TT!.jl")
    global SparseStaticResid! =
        load_dynare_function("$(function_root)SparseStaticResid!.jl")
    global SparseStaticResidTT! =
        load_dynare_function("$(function_root)SparseStaticResidTT!.jl")
    global SparseDynamicParametersDerivatives! =
        load_dynare_function("$(function_root)DynamicParamsDerivs.jl")
    global SparseStaticParametersDerivatives! =
        load_dynare_function("$(function_root)StaticParamsDerivs.jl")
    global steady_state!
    (steady_state!, analytical_variables) =
        load_steady_state_function("$(function_root)SteadyState2.jl")
    global dynamic_auxiliary_variables! =
        load_dynare_function("$(function_root)DynamicSetAuxiliarySeries.jl")
    global static_auxiliary_variables! =
        load_dynare_function("$(function_root)SetAuxiliaryVariables.jl")
    return analytical_variables
end


nearbyint(x::T) where {T<:Real} =
    (abs((x) - floor(x)) < abs((x) - ceil(x)) ? floor(x) : ceil(x))

function get_power_deriv(x::T, p::T, k::Int64) where {T<:Real}
    if (abs(x) < 1e-12 && p > 0 && k > p && abs(p - nearbyint(p)) < 1e-12)
        return 0.0
    else
        dxp = x^(p - k)
        for i = 1:k
            dxp *= p
            p -= 1
        end
        return dxp
    end
end

function dynamic_resid!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicResid!(T, residual, y, x, params, steady_state)
    return nothing
end

dynamic_derivatives!(
    T::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real}
) = dynamic_derivatives!(T, g1.nzval, y, x, params, steady_state)

function dynamic_derivatives!(
    T::AbstractVector{<:Real},
    nzval::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG1!(T, nzval, y, x, params, steady_state)
    return nothing
end

function dynamic_derivatives2!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG2TT!(T, y, x, params, steady_state)
    SparseDynamicG2!(T, g2.nzval, y, x, params, steady_state)
    return nothing
end

function dynamic_derivatives3!(
    T::Vector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    g3::AbstractMatrix{<:Real},
    y::Vector{<:Real},
    x::AbstractVector{<:Real},
    params::Vector{<:Real},
    steady_state::Vector{<:Real},
)
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG2TT!(T, y, x, params, steady_state)
    SparseDynamicG3TT!(T, y, x, params, steady_state)
    SparseDynamicG3!(T, g3.nzval, y, x, params, steady_state)
    return nothing
end

function dynamic!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    dynamic_resid!(T, residual, y, x, params, steady_state)
    return nothing
end

function dynamic!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    dynamic_resid!(T, residual, y, x, params, steady_state)
    dynamic_derivatives!(T, g1, y, x, params, steady_state)
    return nothing
end

function dynamic!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    dynamic_resid!(T, residual, y, x, params, steady_state)
    dynamic_derivatives!(T, residual, g1, y, x, params, steady_state)
    dynamic_derivatives2!(T, residual, g1, g2, y, x, params, steady_state)
    return nothing
end

function dynamic!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    g3::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
    steady_state::AbstractVector{<:Real},
)
    dynamic_resid!(T, residual, y, x, params, steady_state)
    dynamic_derivatives!(T, residual, g1, y, x, params, steady_state)
    dynamic_derivatives2!(T, residual, g1, g2, y, x, params, steady_state)
    dynamic_derivatives3!(T, residual, g1, g2, g3, y, x, params, steady_state)
    return nothing
end

function static_resid!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    SparseStaticResidTT!(T, y, x, params)
    SparseStaticResid!(T, residual, y, x, params)
    return nothing
end

function static_derivatives!(
    T::Vector{<:Real},
    g1::AbstractMatrix{<:Real},
    y::Vector{<:Real},
    x::AbstractVector{<:Real},
    params::Vector{<:Real},
)
    SparseStaticResidTT!(T, y, x, params)
    SparseStaticG1TT!(T, y, x, params)
    SparseStaticG1!(T, g1.nzval, y, x, params)
    return nothing
end

function static_derivatives2!(
    T::Vector{<:Real},
    g2::AbstractMatrix{<:Real},
    y::Vector{<:Real},
    x::AbstractVector{<:Real},
    params::Vector{<:Real},
)
    SparseStaticResidTT!(T, y, x, params)
    SparseStaticG1TT!(T, y, x, params)
    SparseStaticG2TT!(T, y, x, params)
    SparseStaticG2!(T, g2.nzval, y, x, params)
    return nothing
end

function static_derivatives3!(
    T::Vector{<:Real},
    g3::AbstractMatrix{<:Real},
    y::Vector{<:Real},
    x::AbstractVector{<:Real},
    params::Vector{<:Real},
)
    SparseStaticResidTT!(T, y, x, params)
    SparseStaticG1TT!(T, y, x, params)
    SparseStaticG2TT!(T, y, x, params)
    SparseStaticG3TT!(T, y, x, params)
    SparseStaticG3!(T, g3.nzval, y, x, params)
    return nothing
end

function static!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    static_resid!(T, residual, y, x, params)
    return nothing
end

function static!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    static_resid!(T, residual, y, x, params)
    static_derivatives!(T, g1, y, x, params)
    return nothing
end

function static!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    static_resid!(T, residual, y, x, params)
    static_derivatives!(T, g1, y, x, params)
    static_derivatives2!(T, g2, y, x, params)
    return nothing
end

function static!(
    T::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    g1::AbstractMatrix{<:Real},
    g2::AbstractMatrix{<:Real},
    g3::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    static_resid!(T, residual, y, x, params)
    static_derivatives!(T, g1, y, x, params)
    static_derivatives2!(T, g2, y, x, params)
    static_derivatives3!(T, g3, y, x, params)
    return nothing
end


function load_dynare_function(modname::String; head = 1, tail = 0)::Function
    if isfile(modname)
        fun = readlines(modname)
        return (@RuntimeGeneratedFunction(Meta.parse(join(fun[head:end-tail], "\n"))))
    else
        return (x...) -> nothing
    end

end

function load_set_dynamic_auxiliary_variables(modelname::String)
    source = []
    functionstart = false
    for line in readlines("$(modelname)DynamicSetAuxiliarySeries.jl", keep = true)
        if startswith(line, "function")
            functionstart = true
        end
        if functionstart
            push!(source, line)
            if startswith(line, "end")
                push!(source, "end")
                break
            end
        end
    end
    exp1 = Meta.parse(join(source, "\n"))
    return (@RuntimeGeneratedFunction(exp1))
end

function load_steady_state_function(modname::String)
    if isfile(modname)
        fun = readlines(modname)
        expr = Meta.parse(join(fun, "\n"))
        analytical_variables = get_analytical_variables(expr)
        return (@RuntimeGeneratedFunction(expr), analytical_variables)
    else
        return (nothing, [])
    end
end

function get_analytical_variables(expr::Expr)
    block = expr.args[2].args[3].args[3]
    @assert  block.head == :block
    indices = Int64[]
    for a in block.args
        if (isa(a, Expr)
            && a.head == :(=)
            && isa(a.args[1], Expr)
            && a.args[1].args[1] == :ys_)
            push!(indices, a.args[1].args[2])
        end
    end

    return sort(unique(indices))
end

end # end module
