using MyFirstPackage
using Documenter

DocMeta.setdocmeta!(MyFirstPackage, :DocTestSetup, :(using MyFirstPackage); recursive=true)

makedocs(;
    modules=[MyFirstPackage],
    authors="QingyunQian",
    sitename="MyFirstPackage.jl",
    format=Documenter.HTML(;
        canonical="https://QingyunQian.github.io/MyFirstPackage.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/QingyunQian/MyFirstPackage.jl",
    devbranch="main",
)
