using CDGRN
using Documenter

DocMeta.setdocmeta!(CDGRN, :DocTestSetup, :(using CDGRN); recursive=true)

makedocs(;
    modules=[CDGRN],
    authors="Yueh-Hua Tu",
    repo="https://github.com/yuehhua/CDGRN.jl/blob/{commit}{path}#{line}",
    sitename="CDGRN.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuehhua.github.io/CDGRN.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yuehhua/CDGRN.jl",
    devbranch="main",
)
