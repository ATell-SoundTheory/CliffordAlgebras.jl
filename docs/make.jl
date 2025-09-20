push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using Documenter, CliffordAlgebras

makedocs(
    sitename = "CliffordAlgebras.jl Documentation",
    modules = [CliffordAlgebras],
    doctest = true,
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md", 
        "Mathematical Background" => "mathematical_background.md",
        "Performance Guide" => "performance.md",
        "Developer Guide" => "developer.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/ATell-SoundTheory/CliffordAlgebras.jl.git",
    versions = nothing,
)
