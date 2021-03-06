using georectification
using Documenter

DocMeta.setdocmeta!(georectification, :DocTestSetup, :(using georectification); recursive=true)

makedocs(;
    modules=[georectification],
    authors="john.louis.waczak@gmail.com",
    repo="https://github.com/mi3nts/georectification.jl/blob/{commit}{path}#{line}",
    sitename="georectification.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mi3nts.github.io/georectification.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mi3nts/georectification.jl",
    devbranch="main",
)
