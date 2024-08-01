using PATHSolver, Dynare
using Test

context = @dynare "models/global/irbc" "-DN=2" "-DMETHOD=\"NewtonRaphson()\"";
context = @dynare "models/global/irbc" "-DN=2" "-DMETHOD=\"TrustRegion()\"";
context = @dynare "models/global/irbc_irr" "-DN=2" "-DSOLVER=\"NLsolver\"";
context = @dynare "models/global/irbc_irr" "-DN=2" "-DSOLVER=\"PATHSolver\"";
nothing



