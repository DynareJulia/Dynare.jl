using Documenter, Dynare

makedocs(
    sitename="Dynare.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    # doctest = false,
    pages=[
        "Home" => "index.md",
        "Installation and Configuration" => "installation-and-configuration.md",
        "Running Dynare" => "running-dynare.md",
        "The Model File" => "the-model-file.md",
    ],
)

#if get(ENV, "CI", nothing) == "true"
#    deploydocs(; repo="github.com/DynareJulia/Dynare.jl.git", push_preview=true)
#end
