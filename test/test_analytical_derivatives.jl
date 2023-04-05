using Dynare

context = @dynare "models/analytical_derivatives/fs2000_sa.mod" "params_derivs_order=1" "notmpterms";

params = context.work.params
steady_state = context.results.model_results[1].trends.endogenous_steady_state
exogenous = context.results.model_results[1].trends.exogenous_steady_state

(rp, gp) = Dynare.DFunctions.SparseStaticParametersDerivatives!(steady_state, exogenous, params)
