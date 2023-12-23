using Documenter, Dynare

makedocs(
    sitename="Dynare.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # doctest = false,
    pages=[
        "Home" => "index.md",
        "Installation and Configuration" => "installation-and-configuration.md",
        "Running Dynare" => "running-dynare.md",
        "Model File" => [
            "Syntax elements" => "model-file/syntax-elements.md",
            "Variables and parameters declaration" => "model-file/variable-declarations.md",
            "Model declaration" => "model-file/model-declaration.md",
#            "Initial and terminal conditions" => "model-file/initial-terminal-conditions.md",
            "Steady state" => "model-file/steady-state.md",
            "Shocks on exgogenous variables" => "model-file/shocks.md",
            "Deterministic simulations" => "model-file/deterministic-simulations.md",
            "Local approximation" => "model-file/local-approxiation.md",
            "State space, filtering and smoothing" => "model-file/filtersmoother.md",
            "Estimation" => "model-file/estimation.md",
            "Forecasting" => "model-file/forecasting.md",
#            "Optimal policy" => "model-file/optimal-policy.md",
        ],
        "Macroprocessing language" => "macroprocessor.md"     
    ],
    pagesonly = true,
    
)

deploydocs(
    repo="github.com/DynareJulia/Dynare.jl.git", push_preview=true
)
