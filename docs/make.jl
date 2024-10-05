using Documenter, DocumenterCitations, Dynare

# A flag to check if we are running in a GitHub action.
const _IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Pass --pdf to build the PDF. On GitHub actions, we always build the PDF.
const _PDF = findfirst(isequal("--pdf"), ARGS) !== nothing || _IS_GITHUB_ACTIONS
_PAGES =[
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
        "Reporting"=> "model-file/reporting.md",
        #            "Optimal policy" => "model-file/optimal-policy.md",
    ],
    "Macroprocessing language" => "macroprocessor.md",
    "References" => "references.md",
]

# Needed to make Documenter think that there is a PDF in the right place when
# link checking. Inn production we replace this by running the LaTeX build.
write(joinpath(@__DIR__, "src", "Dynare.pdf"), "")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
makedocs(
    sitename="Dynare.jl",
    format=Documenter.HTML(;
                           prettyurls = get(ENV, "CI", nothing) == "true",
                           assets = String["assets/citations.css"]),
    # doctest = false,
    pages = _PAGES,
    pagesonly = true,
    plugins = [bib],
    
)

latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"
bib1 = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
makedocs(
    sitename = "Dynare",
    format = Documenter.LaTeX(; platform = latex_platform),
    build = "latex_build",
    pages = _PAGES,
    pagesonly = true,
    debug = true,
    plugins = [bib1]
    )
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
     cp(
        joinpath(@__DIR__, "latex_build", "Dynare.pdf"),
        joinpath(@__DIR__, "build", "Dynare.pdf");
        force = true,
)

deploydocs(
    repo="github.com/DynareJulia/Dynare.jl.git", push_preview=true
)
