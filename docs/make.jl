using CDGRNs
using Documenter

DocMeta.setdocmeta!(CDGRNs, :DocTestSetup, :(using CDGRNs); recursive=true)

makedocs(;
    modules=[CDGRNs],
    authors="Yueh-Hua Tu",
    repo="https://github.com/yuehhua/CDGRNs.jl/blob/{commit}{path}#{line}",
    sitename="CDGRNs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuehhua.github.io/CDGRNs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yuehhua/CDGRNs.jl",
    devbranch="main",
)
