module DFunctions

using RuntimeGeneratedFunctions
using SparseArrays
using StatsFuns
using TimeDataFrames

RuntimeGeneratedFunctions.init(@__MODULE__)

function load_model_functions(modelname::String)
    function_root = "$(modelname)/model/julia/"
    @show function_root
    @show "$(function_root)SparseDynamicG1!.jl"
    global SparseDynamicG1! = load_dynare_function("$(function_root)SparseDynamicG1!.jl") 
    global SparseDynamicG1TT! = load_dynare_function("$(function_root)SparseDynamicG1TT!.jl")  
    global SparseDynamicG2! = load_dynare_function("$(function_root)SparseDynamicG2!.jl")  
    global SparseDynamicG2TT! = load_dynare_function("$(function_root)SparseDynamicG2TT!.jl")  
    global SparseDynamicG3! = load_dynare_function("$(function_root)SparseDynamicG3!.jl")  
    global SparseDynamicG3TT! = load_dynare_function("$(function_root)SparseDynamicG3TT!.jl")  
    global SparseDynamicResid! = load_dynare_function("$(function_root)SparseDynamicResid!.jl") 
    global SparseDynamicResidTT! = load_dynare_function("$(function_root)SparseDynamicResidTT!.jl") 
    global SparseStaticG1! = load_dynare_function("$(function_root)SparseStaticG1!.jl") 
    global SparseStaticG1TT! = load_dynare_function("$(function_root)SparseStaticG1TT!.jl") 
    global SparseStaticG2! = load_dynare_function("$(function_root)SparseStaticG2!.jl") 
    global SparseStaticG2TT! = load_dynare_function("$(function_root)SparseStaticG2TT!.jl") 
    global SparseStaticG3! = load_dynare_function("$(function_root)SparseStaticG3!.jl") 
    global SparseStaticG3TT! = load_dynare_function("$(function_root)SparseStaticG3TT!.jl") 
    global SparseStaticResid! = load_dynare_function("$(function_root)SparseStaticResid!.jl")  
    global SparseStaticResidTT! = load_dynare_function("$(function_root)SparseStaticResidTT!.jl")
    global SparseDynamicParametersDerivatives! = load_dynare_function("$(modelname)DynamicParamsDerivs.jl", head = 8, tail = 1)
    global SparseStaticParametersDerivatives! = load_dynare_function("$(modelname)StaticParamsDerivs.jl", head = 8, tail = 1)
    global steady_state! = load_dynare_function("$(modelname)SteadyState2.jl", head = 8, tail = 1)
    global dynamic_auxiliary_variables! = load_dynare_function("$(modelname)DynamicSetAuxiliarySeries.jl", head = 3, tail =1 )
    global static_auxiliary_variables! = load_dynare_function("$(modelname)SetAuxiliaryVariables.jl", head = 3, tail =1 )
    return nothing
end
    
    
nearbyint(x::T) where T <: Real  = (abs((x)-floor(x)) < abs((x)-ceil(x)) ? floor(x) : ceil(x))

function get_power_deriv(x::T, p::T, k::Int64) where T <: Real
    if (abs(x) < 1e-12 && p > 0 && k > p && abs(p-nearbyint(p)) < 1e-12 )
        return 0.0
    else
        dxp = x^(p-k)
        for i = 1:k
	     dxp *= p
	     p -= 1
	 end
	 return dxp
    end
end

function dynamic!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractVector{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real})
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicResid!(T, residual, y, x, params, steady_state)
    return nothing
end

function dynamic!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicResid!(T, residual, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG1!(T, g1.nzval, y, x, params, steady_state)
    return nothing
end

function dynamic!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicResid!(T, residual, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG1!(T, g1.nzval, y, x, params, steady_state)
    SparseDynamicG2TT!(T, y, x, params, steady_state)
    SparseDynamicG2!(T, g2.nzval, y, x, params, steady_state)
    return nothing
end

function dynamic!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real}, g3::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseDynamicResidTT!(T, y, x, params, steady_state)
    SparseDynamicResid!(T, residual, y, x, params, steady_state)
    SparseDynamicG1TT!(T, y, x, params, steady_state)
    SparseDynamicG1!(T, g1.nzval, y, x, params, steady_state)
    SparseDynamicG2TT!(T, y, x, params, steady_state)
    SparseDynamicG2!(T, g2.nzval, y, x, params, steady_state)
    SparseDynamicG3TT!(T, y, x, params, steady_state)
    SparseDynamicG3!(T, g3.nzval, y, x, params, steady_state)
    return nothing
end

function static!(T::Vector{<: Real}, residual::AbstractVector{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseStaticResidTT!(T, y, x, params, steady_state)
    SparseStaticResid!(T, residual, y, x, params, steady_state)
    return nothing
end

function static!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseStaticResidTT!(T, y, x, params, steady_state)
    SparseStaticResid!(T, residual, y, x, params, steady_state)
    SparseStaticG1TT!(T, y, x, params, steady_state)
    SparseStaticG1!(T, g1.nzval, y, x, params, steady_state)
    return nothing
end

function static!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseStaticResidTT!(T, y, x, params, steady_state)
    SparseStaticResid!(T, residual, y, x, params, steady_state)
    SparseStaticG1TT!(T, y, x, params, steady_state)
    SparseStaticG1!(T, g1,nzval, y, x, params, steady_state)
    SparseStaticG2TT!(T, y, x, params, steady_state)
    SparseStaticG2!(T, g2.nzval, y, x, params, steady_state)
    return nothing
end

function static!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real}, g3::AbstractMatrix{<: Real},
                  y::Vector{<: Real}, x::AbstractVector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseStaticResidTT!(T, y, x, params, steady_state)
    SparseStaticResid!(T, residual, y, x, params, steady_state)
    SparseStaticG1TT!(T, y, x, params, steady_state)
    SparseStaticG1!(T, g1.nzval, y, x, params, steady_state)
    SparseStaticG2TT!(T, y, x, params, steady_state)
    SparseStaticG2!(T, g2.nzval, y, x, params, steady_state)
    SparseStaticG3TT!(T, y, x, params, steady_state)
    SparseStaticG3!(T, g3.nzval, y, x, params, steady_state)
    return nothing
end

function load_dynare_function(modname::String; head=1, tail=0)::Function
    if isfile(modname)
        fun = readlines(modname)
        return (@RuntimeGeneratedFunction(Meta.parse(join(fun[head:end-tail], "\n"))))
    else
        return (x...) -> nothing
    end
            
end


end # end module
