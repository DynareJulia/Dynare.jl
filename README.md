# WORK IN PROGRESS

## Installation 

```
using Pkg
pkg"add Dynare"
```
## Update
If you have already a version of Dynare installed (not in development mode):
using Pkg
pkg"update Dynare"

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
