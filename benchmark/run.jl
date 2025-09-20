#!/usr/bin/env julia
# Lightweight runner for the benchmark suite.
# Usage: julia --project=benchmark benchmark/run.jl

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Ensure the local package is available in this environment
try
	Base.find_package("CliffordAlgebras") === nothing && error()
catch
	repo_root = normpath(joinpath(@__DIR__, ".."))
	Pkg.develop(PackageSpec(path = repo_root))
end

using BenchmarkTools
include(joinpath(@__DIR__, "benchmarks.jl"))

# Run a small warmup and then the full suite
println("Running CliffordAlgebras.jl benchmarks...")
warmup(SUITE)
results = run(SUITE; verbose = true)
show(results)
println()
