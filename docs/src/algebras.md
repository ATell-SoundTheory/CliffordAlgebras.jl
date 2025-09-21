# Algebra Catalog

This page lists predefined algebra aliases available in `CliffordAlgebras.jl`, along with their signatures and quick facts.

| Alias(es) | Signature (p,q,r) | Base symbols | Notes |
|---|---:|---|---|
| `:Cl2` | (2,0,0) | `(:e1, :e2)` | 2D Euclidean; complex-plane-like bivector `e1e2` |
| `:Cl3` | (3,0,0) | `(:e1, :e2, :e3)` | 3D Euclidean; rotors from bivectors |
| `:Spacetime`, `:STA` | (1,3,0) | `(:t, :x, :y, :z)` | Minkowski spacetime |
| `:Complex`, `:ℂ` | (0,1,0) | `(:i,)` | Complex numbers as a 2D algebra |
| `:Quaternions`, `:ℍ` | (0,2,0) | `(:i, :j)` | Even subalgebra of Cl(3,0,0) isomorphic to quaternions |
| `:Hyperbolic`, `:Hyper` | (1,0,0) | `(:j,)` | Hyperbolic numbers |
| `:Dual`, `:Grassmann` | (0,0,1) | `(:ε,)` | Dual numbers (nilpotent) |
| `:Grassmann2D`, `:G2` | (0,0,2) | `(:ε₁, :ε₂)` | 2D Grassmann (both square to 0) |
| `:Grassmann3D`, `:G3` | (0,0,3) | `(:ε₁, :ε₂, :ε₃)` | 3D Grassmann |
| `:PGA2D`, `:Projective2D`, `:Plane2D` | (2,0,1) | `(:e1, :e2, :e0)` | Projective GA (2D); `e0^2=0` |
| `:PGA3D`, `:Projective3D`, `:Plane3D` | (3,0,1) | `(:e1, :e2, :e3, :e0)` | Projective GA (3D); `e0` is null |
| `:CGA2D`, `:Conformal2D` | (3,1,0) | `(:e1, :e2, :e₊, :e₋)` | Conformal GA (2D) |
| `:CGA3D`, `:Conformal3D` | (4,1,0) | `(:e1, :e2, :e3, :e₊, :e₋)` | Conformal GA (3D) |
| `:DCGA3D`, `:DoubleConformal3D` | (6,2,0) | — | Double conformal (3D) |
| `:TCGA3D`, `:TripleConformal3D` | (9,3,0) | — | Triple conformal (3D) |
| `:DCGSTA`, `:DoubleConformalSpacetime` | (4,8,0) | `(:t₁, :t₂, :e₊₁, :e₊₂, :x₁, :x₂, :y₁, :y₂, :z₁, :z₂, :e₋₁, :e₋₂)` | DCG for spacetime |
| `:QCGA`, `:QuadricConformal` | (9,6,0) | — | Quadric conformal |

Tip: Use `signaturetable(stdout, algebra)` to view per-basis signatures and `cayleytable(stdout, algebra)` for the full multiplication table.

## Quick facts and size

For an algebra with signature (p,q,r) and order n = p+q+r:

- Dimension of the full algebra: 2^n elements.
- Number of k-vectors (grade k): C(n,k).
- Pseudoscalar I has grade n and I^2 = character(algebra) ∈ {+1,-1,0}.
- Null basis elements (r > 0) square to 0 and model points at infinity or conformal components.

Examples:
- Cl(3): n=3, dim=8, grades per k: 1, 3, 3, 1.
- PGA3D (3,0,1): n=4, dim=16, with one null basis e0.
- CGA3D (4,1,0): n=5, dim=32, two lightlike directions from e₊, e₋.

## Typical use cases

- Cl(2), Cl(3): Euclidean plane/space; basic rotations, rigid body kinematics (rotors from bivectors).
- STA (Cl(1,3,0)): Relativistic spacetime computations.
- PGA2D/PGA3D: Projective geometry for graphics/robotics; lines/planes at infinity via null basis `e0`.
- CGA2D/CGA3D: Conformal geometry for points, circles/spheres, and conformal transforms.
- ℂ / ℍ: Complex/quaternion arithmetic embedded in geometric algebra contexts.

## Examples

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> io = IOBuffer(); signaturetable(io, pga); true
true

julia> e0 = basevector(pga, :e0);  # null basis vector

julia> scalar(e0*e0)
0

julia> cga = CliffordAlgebra(:CGA3D);

julia> eplus = basevector(cga, :e₊); eminus = basevector(cga, :e₋);

julia> (scalar(eplus*eplus), scalar(eminus*eminus))
(1, -1)
```
