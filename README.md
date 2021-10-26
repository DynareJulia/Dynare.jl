# WORK IN PROGRESS

## Installation 

```
using Pkg
pkg"registry add https://git.dynare.org/DynareJulia/DynareRegistry.git"
pkg"add Dynare"
```
in case of difficulty with the installation of ``Plots``, try
```
pkg"build Plots"
pkg"build Dynare"
```
## Running Dynare


Example (to be run in the directory ``./Dynare.jl``):
```
using Dynare
context = @dynare "./test/models/example1/example1.mod";
```
The results are in the ``context`` structure.
