using MomentFittingQuadrature
using Documenter

DocMeta.setdocmeta!(MomentFittingQuadrature, :DocTestSetup, :(using MomentFittingQuadrature); recursive=true)

makedocs(;
    modules=[MomentFittingQuadrature],
    authors="Michał Mika and contributors",
    sitename="MomentFittingQuadrature.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/MomentFittingQuadrature.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/MomentFittingQuadrature.jl",
    devbranch="main",
)
