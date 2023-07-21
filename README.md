# WORK IN PROGRESS

## Requirements
Julia 64bit >= v1.6

## Installation 

```
using Pkg
pkg"add Dynare"
```
## Update
If you already have a version of Dynare installed (not in development mode):
```
using Pkg
pkg"update"
```
## Running Dynare

Example (to be run in the directory ``Dynare.jl``):
```
using Dynare
context = @dynare "./test/models/example1/example1.mod";
```
The results are in the ``context`` structure.

## Supported Dynare instructions

1. calib_smoother
1. check
1. deterministic_trends
1. endval
1. histval
1. initval
1. perfect\_foresight\_setup (only some options)
1. perfect\_foresight\_solver (only some options, includind lmmcp)
1. planner_objective
1. ramsey\_model
1. shocks
1. steady (including numerical solution)
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

## Extensions
1. PATH algorithm for MCP problems
  - You need to load PATHSolver before running Dynare (once per Julia
    session)
```
	using PATHSolver
    context = @dynare ...
```
1. PARDISO is high-performance solver for very large sparse linear
   systems. It can be used with perfect foresight simulation of very
   large models. You need to load MKL and Pardiso before running
   Dynare (once per Julia session)
```
   using MKL, Pardiso
   context = @dynare ...
```

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

### Under Windows
1. Create a new folder in your ``.julia`` directory called ``config``
1. In the ``.julia\config`` folder create a new text file: ``Right click  → New → Text Document``
1. Inside the new text document, type:

> ENV["PATH_LICENSE_STRING"] = "*licence number provided in the above
> link*"

1. In the text document, go ``File → Save as``
1. In the dropdown ``Save as type`` option at the bottom of the pop-up window select ``All files``
1. In the field ``File name`` write ``startup.jl``
1. Press ``Save``

### In VScode
1. Go to ``File → Open Folder``, navigate through your folders and choose the one where you want to create the ``startup.jl`` file
1. Go to ``File → New Text File``
1. In the new text file, type
> ENV["PATH_LICENSE_STRING"] = "*licence number provided in the above
> link*"

1. Go to ``File → Save``
1. In the dropdown ``Save as type`` option at the bottom of the pop-up window select ``All files``
1. In the field ``File name`` write ``startup.jl``.
1. Press ``Save``
