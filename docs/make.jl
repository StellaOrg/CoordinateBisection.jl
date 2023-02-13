using CoordinateBisection
using Documenter

DocMeta.setdocmeta!(CoordinateBisection, :DocTestSetup, :(using CoordinateBisection); recursive=true)

makedocs(;
    modules=[CoordinateBisection],
    authors="Andrei-Leonard Nicusan <a.l.nicusan@gmail.com> and contributors",
    repo="https://github.com/anicusan/CoordinateBisection.jl",
    sitename="CoordinateBisection.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://anicusan.github.io/CoordinateBisection.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anicusan/CoordinateBisection.jl",
    devbranch="main",
)
