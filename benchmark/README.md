# Benchmarks

This folder contains a small benchmark suite for CliffordAlgebras.jl.

Two ways to run:

- Quick run (local):
  - julia --project=benchmark benchmark/run.jl

- With PkgBenchmark (compare against a baseline):
  - julia --project=benchmark -e 'using PkgBenchmark, CliffordAlgebras; benchmarkpkg(CliffordAlgebras)'

Scenarios included:
- Cl(2): geometric product, wedge, rotor sandwich
- Cl(3): geometric product, dot product, reverse, matrix conversion
- Internals: multiplication table cache retrieval

Notes:
- First run may include precompilation; subsequent runs are more stable.
- For fair comparisons, close other apps and pin Julia threads if needed (`JULIA_NUM_THREADS`).
