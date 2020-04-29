# WORK IN PROGRESS

## Installation 

### Adding packages 
```
using Pkg
Pkg.add(PackageSpec(url="https://git.dynare.org/julia-packages/fastlapackinterface.jl.git"))
Pkg.add(PackageSpec(url="https://git.dynare.org/julia-packages/polynomialmatrixequations.jl.git"))
Pkg.add(PackageSpec(url="https://git.dynare.org/julia-packages/linearrationalexpectations.jl.git"))
Pkg.add(PackageSpec(url="https://git.dynare.org/julia-packages/dynare.jl.git"))
```
or, interactively,
```
]add https://git.dynare.org/julia-packages/fastlapackinterface.jl.git
]add https://git.dynare.org/julia-packages/polynomialmatrixequations.jl.git
]add https://git.dynare.org/julia-packages/linearrationalexpectations.jl.git
]add https://git.dynare.org/julia-packages/dynare.jl.git
```

Until a proper ``build.jl`` is added, one must confiure DYNARE_ROOT in ``./src/parser/DynarePreporcessor.jl`` to point to the installation of a current unstable version of Dynare or Dynare preprocessor
