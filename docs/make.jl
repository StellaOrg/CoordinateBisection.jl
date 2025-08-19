using CoordinateBisection
using Documenter


makedocs(;
    modules=[CoordinateBisection],
    authors="Andrei-Leonard Nicusan",
    sitename="CoordinateBisection.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/StellaOrg/CoordinateBisection.jl",     # TODO: update this to docs webpage
        edit_link="main",
        assets=String[],
        sidebar_sitename=false,

        # Only create web pretty-URLs on the CI
        prettyurls=get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Overview" => "index.md",
        # "Benchmarks" => "benchmarks.md",
        # "Performance Tips" => "performance.md",
        # "Manual" =>[
        #     "Using Different Backends" => "api/using_backends.md",
        #     "General Loops" => "api/foreachindex.md",
        #     "Map" => "api/map.md",
        #     "Sorting" => "api/sort.md",
        #     "Reduce" => "api/reduce.md",
        #     "MapReduce" => "api/mapreduce.md",
        #     "Accumulate" => "api/accumulate.md",
        #     "Binary Search" => "api/binarysearch.md",
        #     "Predicates" => "api/predicates.md",
        #     "Arithmetics" => "api/arithmetics.md",
        #     "Custom Structs" => "api/custom_structs.md",
        #     "Task Partitioning" => "api/task_partition.md",
        #     "Utilities" => "api/utilities.md",
        # ],
        # "Testing" => "testing.md",
        # "Debugging Kernels" => "debugging.md",
        # "Roadmap" => "roadmap.md",
        # "References" => "references.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="https://github.com/StellaOrg/CoordinateBisection.jl",
    devbranch="main",
)
