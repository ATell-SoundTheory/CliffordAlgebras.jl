using Documenter, CliffordAlgebras

makedocs(
    sitename = "CliffordAlgebras.jl Documentation",
    modules = [CliffordAlgebras]
)

deploydocs(
    repo = "github.com/ATell-SoundTheory/CliffordAlgebras.jl.git",
    versions = nothing,
)
