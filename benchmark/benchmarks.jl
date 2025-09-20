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

SUITE["Cl2"]["fat_dot"] = @benchmarkable begin
    a ⋅ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["scalar_product"] = @benchmarkable begin
    a ⋆ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["left_contraction"] = @benchmarkable begin
    a ⨼ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["right_contraction"] = @benchmarkable begin
    a ⨽ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["commutator"] = @benchmarkable begin
    a ×₋ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["anticommutator"] = @benchmarkable begin
    a ×₊ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 2.0*e1 + 3.0*e2; b = 4.0*e1 + 5.0*e2)

SUITE["Cl2"]["regressive_product"] = @benchmarkable begin
    a ∨ b
end setup=(cl2 = CliffordAlgebra(2); e1 = cl2.e1; e2 = cl2.e2; a = 1 + e1; b = 1 + e2)

SUITE["Cl2"]["dual_norm"] = @benchmarkable begin
    norm(dual(mv))
end setup=(cl2 = CliffordAlgebra(2); mv = 1.0 + 2.0*cl2.e1 + 3.0*cl2.e1e2)

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

SUITE["Cl3"]["bivector_wedge_vector"] = @benchmarkable begin
    (e1 ∧ e2) ∧ e3
end setup=(cl3 = CliffordAlgebra(3); e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3)

SUITE["Cl3"]["contractions"] = @benchmarkable begin
    mv1 ⨼ mv2 + mv1 ⨽ mv2
end setup=(cl3 = CliffordAlgebra(3); e1 = cl3.e1; e2 = cl3.e2; e3 = cl3.e3; mv1 = e1 + 2e2 + 3e3; mv2 = 2e1 + e2 - e3)

SUITE["Cl3"]["dual_and_grade"] = @benchmarkable begin
    Λᵏ(dual(mv), 1)
end setup=(cl3 = CliffordAlgebra(3); mv = 1.0 + 2.0*cl3.e1 + 3.0*cl3.e2 + 0.5*cl3.e1e2)

# -----------------------------------
# Internals: multiplication table cache
# -----------------------------------
# Measure the (cached) multiplication table retrieval to ensure the cache stays fast.
SUITE["Internals"] = BenchmarkGroup()
SUITE["Internals"]["multtable_cached_Cl4"] = @benchmarkable begin
    CliffordAlgebras.multtable(CA)
end setup=(CA = typeof(CliffordAlgebra(4)); _ = CliffordAlgebras.multtable(CA))

SUITE["Internals"]["multtable_cached_Cl5"] = @benchmarkable begin
    CliffordAlgebras.multtable(CA)
end setup=(CA = typeof(CliffordAlgebra(5)); _ = CliffordAlgebras.multtable(CA))

# --------------------------
# Larger algebras: Cl(5)
# --------------------------
SUITE["Cl5"] = BenchmarkGroup()
SUITE["Cl5"]["gp_vectors"] = @benchmarkable begin
    v1 * v2
end setup=(cl5 = CliffordAlgebra(5); v1 = 1*cl5.e1 + 2*cl5.e2 + 3*cl5.e3 + 4*cl5.e4 + 5*cl5.e5; v2 = 5*cl5.e1 - 4*cl5.e2 + 3*cl5.e3 - 2*cl5.e4 + 1*cl5.e5)

SUITE["Cl5"]["gp_general"] = @benchmarkable begin
    a * b
end setup=(cl5 = CliffordAlgebra(5); a = 1 + 0.5*cl5.e1 + 0.3*cl5.e2 + 0.2*cl5.e3 + 0.1*cl5.e4 + 0.05*cl5.e5 + 0.7*cl5.e1e2 + 0.4*cl5.e3e4; b = 0.2 + 0.6*cl5.e2 + 0.9*cl5.e5 + 0.3*cl5.e2e3)

SUITE["Cl5"]["reverse_dual_norm"] = @benchmarkable begin
    norm(dual(~a))
end setup=(cl5 = CliffordAlgebra(5); a = 1 + 2*cl5.e1 + 3*cl5.e2 + 4*cl5.e3 + 5*cl5.e4 + 6*cl5.e5)

# --------------------------------
# Larger algebra: PGA3D (Cl(3,0,1))
# --------------------------------
SUITE["PGA3D"] = BenchmarkGroup()
SUITE["PGA3D"]["gp_general"] = @benchmarkable begin
    a * b
end setup=(pga = CliffordAlgebra(:PGA3D); a = 1 + pga.e1 + 2*pga.e2 + 3*pga.e3 + pga.e1e2; b = 2 + 3*pga.e2 + pga.e3 + 2*pga.e1e3)

SUITE["PGA3D"]["sandwich"] = @benchmarkable begin
    R ≀ v
end setup=(pga = CliffordAlgebra(:PGA3D); v = pga.e1 + 2*pga.e2; B = 0.1*pga.e1e2; R = exp(B))

# -------------------------------
# Larger algebra: CGA3D (Cl(4,1))
# -------------------------------
SUITE["CGA3D"] = BenchmarkGroup()
SUITE["CGA3D"]["gp_vectors"] = @benchmarkable begin
    v1 * v2
end setup=(cga = CliffordAlgebra(:CGA3D); v1 = cga.e1 + 2*cga.e2 + 3*cga.e3 + cga.e4 + cga.e5; v2 = 2*cga.e1 - cga.e2 + 0.5*cga.e3 + cga.e4 - cga.e5)

SUITE["CGA3D"]["bivector_exp"] = @benchmarkable begin
    exp(B)
end setup=(cga = CliffordAlgebra(:CGA3D); B = 0.05*(cga.e1e2 + cga.e2e3 + cga.e3e1))

SUITE["CGA3D"]["matrix_conversion"] = @benchmarkable begin
    matrix(a)
end setup=(cga = CliffordAlgebra(:CGA3D); a = 1 + cga.e1 + 2*cga.e2 + 3*cga.e3 + cga.e1e2 + cga.e2e3)
