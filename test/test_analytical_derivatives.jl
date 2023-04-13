module TestDerivatives
using Dynare #Use analytical_derivatives branch
using LinearAlgebra
using FastLapackInterface
using Test

#Model
context = @dynare "models/analytical_derivatives/fs2000_sa.mod" "params_derivs_order=1" "notmpterms";
model = context.models[1]
results = context.results.model_results[1]
wss = Dynare.StaticWs(context)
params = context.work.params
steadystate = results.trends.endogenous_steady_state
endogenous = repeat(steadystate, 3)
exogenous = results.trends.exogenous_steady_state 

#Get StaticJacobian
df_dx = Dynare.get_static_jacobian!(wss, params, steadystate, exogenous, model) 

#Get StaticJacobianParams
(df_dp, gp) = Dynare.DFunctions.SparseStaticParametersDerivatives!(steadystate, exogenous, params)
 
#Get Jacobian matrix using the inverse of df_dx
function Derivatives(df_dx, df_dp)
    dense_df_dx = Matrix(df_dx)

    if det(dense_df_dx) == 0
        error("The matrix is not invertible")
    end

    #Get dx_dp
    df_dx_inv = dense_df_dx \ I#inv(dense_df_dx)
    return -df_dx_inv * df_dp
end

#Get Jacobian matrix using FastLapackInterface (dense matrix)
function Derivatives2(df_dx, df_dp)
    dense_df_dx = Matrix(df_dx)
    LU_df_dx = LU(LAPACK.getrf!(LUWs(dense_df_dx), dense_df_dx)...)
    if any(diag(LU_df_dx.U) .== 0)
        error("The matrix is not invertible")
    end

    #Get dx_dp
    dx_dp = -(LU_df_dx \ df_dp)
    return dx_dp
end

#Get Jacobian matrix using lu(df_dx) (sparse matrix)
function Derivatives3(df_dx, df_dp)
    LU_df_dx = lu(df_dx)

    if any(diag(LU_df_dx.U) .== 0)
        error("The matrix is not invertible")
    end

    # Get dx_dp
    dx_dp = -(LU_df_dx \ df_dp)
    return dx_dp
end

#Get Jacobian matrix using ldiv! (dense matrix)
function Derivatives4(df_dx, df_dp)
    dense_df_dx = Matrix(df_dx)
    LU_df_dx = lu(dense_df_dx)

    if any(diag(LU_df_dx.U) .== 0)
        error("The matrix is not invertible")
    end

    # Preallocate dx_dp
    n, m = size(df_dp)
    dx_dp = Matrix{Float64}(undef, n, m)

    # Get dx_dp
    for i in 1:m
        dp_col = view(df_dp, :, i)
        sol_col = view(dx_dp, :, i)
        ldiv!(sol_col, LU_df_dx, dp_col)
        sol_col .*= -1
    end

    return dx_dp
end


# @show Derivatives(df_dx, df_dp);
@time sol1 = Derivatives(df_dx, df_dp);
@time sol2 = Derivatives2(df_dx, df_dp);
@time sol3 = Derivatives3(df_dx, df_dp);
@time sol4 = Derivatives4(df_dx, df_dp);
@test sol1 ≈ sol2
@test sol1 ≈ sol3
@test sol1 ≈ sol4 
end # end module   
