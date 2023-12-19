# Installation and configuration

## Software requirements
Dynare is available for all plateforms supported by the current
stable version of Julia
(https://julialang.org/downloads/#supported_platforms). It should also
work with older versions of Julia starting with version 1.6.3

## Installation of Dynare

In Julia, install the Dynare.jl package:
```
using Pkg
pkg"add Dynare"
```

## Using Dynare
In order to start using Dynare in Julia, type
```
using Dynare
```

## Optional extensions

### Pardiso
If you want the solution of very large perfect foresight models and
reduce the memory consumption, use the Pardiso package
(https://github.com/JuliaSparse/Pardiso.jl) and type

```
using MKL, Pardiso
```

### PATHSolver
If youw want to solve perfect foresight models with occasionally
binding constraints use the PATHSolver package
(https://github.com/chkwon/PATHSolver.jl) and type
```
using PATHSolver
```

