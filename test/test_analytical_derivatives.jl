module TestDerivatives
using Dynare #Use analytical_derivatives branch
using LinearAlgebra

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
 
#Get Jacobian matrix
function Derivatives(df_dx, df_dp)
    dense_df_dx = Matrix(df_dx)

    if det(dense_df_dx) == 0
        error("The matrix is not invertible")
    end

    #Get dx_dp
    df_dx_inv = inv(dense_df_dx)
    return -df_dx_inv * df_dp
end

@show Derivatives(df_dx, df_dp)

end # end module   
