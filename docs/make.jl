using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter, CliffordAlgebras

# Ensure doctest blocks have access to exported symbols
DocMeta.setdocmeta!(CliffordAlgebras, :DocTestSetup, :(using CliffordAlgebras); recursive=true)

makedocs(
    sitename = "CliffordAlgebras.jl Documentation",
    modules = [CliffordAlgebras],
    doctest = true,
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md", 
        "Mathematical Background" => "mathematical_background.md",
        "Performance Guide" => "performance.md",
        "Developer Guide" => "developer.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/ATell-SoundTheory/CliffordAlgebras.jl.git",
    devbranch = "master",
    versions = nothing,
)
