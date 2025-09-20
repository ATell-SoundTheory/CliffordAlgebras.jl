"""
Benchmark suite for CliffordAlgebras.jl

Usage with PkgBenchmark (recommended):
  using PkgBenchmark
  benchmarkpkg(CliffordAlgebras)

The suite defines `SUITE` as required by PkgBenchmark.
"""

using BenchmarkTools
using CliffordAlgebras
using LinearAlgebra

const SUITE = BenchmarkTools.BenchmarkGroup()

# ---------------------
# Cl(2): 2D operations
# ---------------------
SUITE["Cl2"] = BenchmarkGroup()
SUITE["Cl2"]["geometric_product"] = @benchmarkable begin
    a * b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["wedge_product"] = @benchmarkable begin
    e1 ∧ e2
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2)

SUITE["Cl2"]["rotor_sandwich"] = @benchmarkable begin
    R ≀ v
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; B = (π/8)*cl2.e1e2; R = exp(B); v = e1)

# ---------------------
# Cl(3): 3D operations
# ---------------------
SUITE["Cl3"] = BenchmarkGroup()
SUITE["Cl3"]["geometric_product"] = @benchmarkable begin
    mv1 * mv2
end setup=(cl3 = CliffordAlgebra(3);
           e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3; e12 = cl3.e1e2;
           mv1 = 1.2 + 2.3*e1 - 0.7*e2 + 3.1*e12 + 0.9*e3;
           mv2 = -0.4 + 1.1*e2 + 2.2*e3 + 0.6*e12)

SUITE["Cl3"]["dot_product"] = @benchmarkable begin
    (e1 + e2) ⋅ (e2 + e3)
end setup=(cl3 = CliffordAlgebra(3); e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3)

SUITE["Cl3"]["reverse"] = @benchmarkable begin
    ~mv
end setup=(cl3 = CliffordAlgebra(3); e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3; e12 = cl3.e1e2; mv = 1.2 + 2.3*e1 - 0.7*e2 + 3.1*e12 + 0.9*e3)

SUITE["Cl3"]["matrix_conversion"] = @benchmarkable begin
    matrix(mv)
end setup=(cl3 = CliffordAlgebra(3); e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3; e12 = cl3.e1e2; mv = 1.2 + 2.3*e1 - 0.7*e2 + 3.1*e12 + 0.9*e3)

# -----------------------------------
# Internals: multiplication table cache
# -----------------------------------
# Measure the (cached) multiplication table retrieval to ensure the cache stays fast.
SUITE["Internals"] = BenchmarkGroup()
SUITE["Internals"]["multtable_cached_Cl4"] = @benchmarkable begin
    CliffordAlgebras.multtable(CA)
end setup=(CA = typeof(CliffordAlgebra(4)); _ = CliffordAlgebras.multtable(CA))
