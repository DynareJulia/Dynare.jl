# WORK IN PROGRESS

## Installation 

```
using Pkg
pkg"add Dynare"
```
## Update
If you have already a version of Dynare installed (not in development mode):
```
using Pkg
pkg"update Dynare"
```
## Running Dynare


Example (to be run in the directory ``./Dynare.jl``):
```
using Dynare
context = @dynare "./test/models/example1/example1.mod";
```
The results are in the ``context`` structure.

## Supported Dynare instructions

1. calib_smoother
1. check
1. deterministic_trends
1. histval
1. initval
1. perfect\_foresight\_setup (only some options)
1. perfect\_foresight\_solver (only some options, includind lmmcp)
1. planner_objective
1. ramsey\_model
1. shocks
1. steady
1. stoch_simul (only order=1)

## Output
1. The ``context`` structure is saved in the directory
   ``<path to modfile>/<modfilenane>/output/<modfilename>.jld2``. It can be loaded with
   ```
   using JLD2
   DD = load("<path to modfile>/<modefilename>/output/<modefilename>.jld2")``
   ```
1. The IRF graphs are saved in ``<path to
   modfile>/<modfilenane>/graphs``
   
## Project web site

- https://www.dynare.org

## PATH licence

Dynare uses the PATH software by S. Dirkse, M.C. Ferris and T. Munson (https://pages.cs.wisc.edu/~ferris/path.html),
as provided by
[PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl) to solve
prefect foresight models with occasionally binding constraints. 

In order to use it, you need to add the free licence available at
https://pages.cs.wisc.edu/~ferris/path/LICENSE in your file
``~/.julia/config/startup.jl``

> ENV["PATH_LICENSE_STRING"] = "*licence number provided in the above
> link*"

