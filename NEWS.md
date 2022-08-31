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
