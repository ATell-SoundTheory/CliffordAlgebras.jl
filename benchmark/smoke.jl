#!/usr/bin/env julia
# Tiny smoke benchmark run intended for CI to catch glaring regressions.
# It runs a very small subset quickly without storing artifacts.

using Pkg
pkgdir = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(@__DIR__)
try
    Pkg.develop(PackageSpec(path = pkgdir))
catch
end
Pkg.instantiate()

using BenchmarkTools
using CliffordAlgebras

cl3 = CliffordAlgebra(3)
e1, e2, e3 = cl3.e1, cl3.e2, cl3.e3

# Minimal set
ops = Dict{String,Function}(
    "gp" => () -> (e1 + 2e2) * (e2 + 3e3),
    "wedge" => () -> (e1 ∧ e2),
    "reverse" => () -> ~(1 + e1 + e2 + e1∧e2),
)

for (name, thunk) in ops
    trial = @benchmarkable thunk() samples=10 evals=1 seconds=0.2
    res = run(trial)
    println(name, ": median = ", BenchmarkTools.median(res))
end
