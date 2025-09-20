# Developer Guide

Welcome! This page is a quick pointer for contributors.

- See `CONTRIBUTING.md` at the repo root for setup, coding style, and PR tips.
- See `CHANGELOG.md` for a summary of notable changes.

## Running the tests

- Unit tests live in `test/runtests.jl`.
- Doctests are enabled and run as part of the docs build.
- Aqua.jl quality checks run from `test/aqua.jl` and are included in `runtests.jl`.

## Useful commands

- Run the package tests locally (uses your default Julia):

```julia
using Pkg
Pkg.activate(temp=true)
Pkg.develop(path="/absolute/path/to/CliffordAlgebras.jl")
Pkg.test("CliffordAlgebras")
```

- Build the docs locally:

```julia
using Pkg
Pkg.activate("docs")
Pkg.instantiate()
include("docs/make.jl")
```

## Notes on precompilation

We include a small `PrecompileTools` workload in `src/precompile.jl` to reduce first-call latency without slowing precompilation. Keep it lightweight; avoid calling heavy or world-age-sensitive paths.

## Reporting issues

Please open an issue with a minimal reproducer and version info (Julia version, OS, package version). Thanks for contributing!
