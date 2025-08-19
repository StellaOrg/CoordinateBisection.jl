using CoordinateBisection
using Documenter


makedocs(;
    modules=[CoordinateBisection],
    authors="Andrei-Leonard Nicusan",
    sitename="CoordinateBisection.jl",
    format=Documenter.HTML(;
        canonical="https://stellaorg.github.io/CoordinateBisection.jl/",
        edit_link="main",
        assets=String[],
        sidebar_sitename=false,

        # Only create web pretty-URLs on the CI
        prettyurls=get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Overview" => "index.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/StellaOrg/CoordinateBisection.jl",
    devbranch="main",
)
