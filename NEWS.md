0.9.3
=====
- add support for purely backward linear models


0.9.2
=====
- fix bug in rwmh_compute!()

0.9.1
=====
- restore printing of Dynare parsing errors

0.9.0
=====
- add normcdf and normpdf
- add periods for calibsmoother and forecast
- add smoothed shocks computation
- fix steady state in filter and smoother
- add forecasting and recursive forecasting
- fix endval
- add posterior mode computation
- add homotopy for steady state 
- add exported functions
- add docstrings
- add documentation

0.8.2
=====
- update DynarePreprocessor_jll to v6.4.0

0.8.1
=====
- conditional simulations with MCP
- add options to perfect_foresight!
- fix to 2nd order computation
- fix to nonstationary models
- fix autocorrelation 

0.8.0
=====
- adding scenario!() and conditional simulations
- adding Julia functions equivalent to some Dynare instructions
  - localapproximation!() -> stoch_simul
  - steadystate!() -> steady
  - perfect_foresight!() -: perfect_foresight_setup and
    perfect_foresight_solver
- adding second order approximation (limited features)
- put PathSolver and PardisoSolver as Julia extensions
- fix occasional problem with static variables and QR decomposition

0.7.5
=====
- fix bug in initialization of perfect foresight simulation

0.7.4
=====
- fix bug in estimation in 0.7.3

0.7.3
=====
- fix initialization with lags and leads on several periods

0.7.2
=====
- fix initialization of perfect foresight
- fix handling of estimated parameters and priors

0.7.1
====
- fix LTS compatibilty
- fix lint errors

0.7.0
=====
- adding estimation module
- replacing TimeDataFrames by AxisArrayTables
- improving simulations storing

0.6.3
=====
- fix compatibility with DynarePreprocessor_jll

0.6.2
=====
- fix initval with parameters
- fix occasionally binding constraints

0.6.1
=====
- fix histval for models with lead or lag exogenous variable
- fix shocks for exogenous variables with steady state different from zero
- add tests for initialization

0.6.0
====
- add support for endval
- fix timing bug in perfect foreisght
- change location of Dynare model functions
- switch to sparse matrices for model derivatives

0.5.7
=====
- add numerical solution for steady state

0.5.6
=====
- fix mcp simulations

0.5.5
=====
- update version of FastLapackInterface, PolynomialMatrixEquations and
  LinearRationalExpectations
  
0.5.4
=====
- add line breaks to model in Latex reports using lstlisting
- provide an additional format where the parameter value is set after
the parameter rather than below

0.5.3
=====
- fix histval for auxiliary variables
- fix histval for prefect foresight models
- fix autocorrelation for stationary models
- don't compute autocorrelation for non stationary models
- add tolf option for steady and perfect foresight
- fix IRFs plotting in several panels (more than 6 variables)
- fix perfect foresight wit no exogenous variable or no shock

0.5.2
=====
- change signature of steady_state_function
- force Int64 numerical constant (necessary for Julia 32bit)
- removing debug instructions

0.5.1
=====
- fix reporting with end of line CR LF

0.5.0
=====
- add plots for TimeDataFrames series
- add accessor functions

0.4.1
=====
- fix ``histval``
- fix ``ramsey_constraints``

0.4.0
=====
- add support for MCP problems in ``perfect_foresight``
- fix preprocessor problem with macOS
- fix backward-looking model simulation
- IRFs are now stored in a vector of TimeDataFrames

0.3.2
=====
- fix graphs in Windows
- fix bugs in initval-file

0.3.1
=====
- update to ExtendedDates v0.1.2
- fix bugs in ``perfect_foresight``

0.3.0
=====
- adding moments and IRFs to stoch_simul
- open automatically image viewer in ``graphs`` directory upon exit
- save ``context`` in a JLD2 file in the ``output`` directory
- fix bug in ``perfect_foresight``

0.2.1
=====
- fixing bug linked to the use of StatsFuns in Dynare generated functions

0.2.0
======
- Moving project from https://git.dynare.org/DynareJulia/Dynare.jl to https://github.com/DynareJulia/Dynare.jl
- Adding project to Julia's General Registry
