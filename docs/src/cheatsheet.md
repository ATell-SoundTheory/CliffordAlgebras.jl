# Cheat Sheet

A quick reference for common tasks, operators, and patterns in CliffordAlgebras.jl.

```@setup cheatsheet
using CliffordAlgebras
using LinearAlgebra
```

```jldoctest cheatsheet
julia> using CliffordAlgebras; using LinearAlgebra; true
true
```

## Algebra creation

```jldoctest cheatsheet
julia> Cl = CliffordAlgebra;  # alias

julia> cl2 = Cl(2); cl3 = Cl(3); sta = Cl(:Spacetime); pga = Cl(:PGA3D); cga = Cl(:CGA3D); true
true
```

## Unicode operators (with ASCII fallbacks)

- Geometric product: `a * b`
- Exterior (wedge): `a ∧ b` — ASCII: `CliffordAlgebras.exteriorprod(a,b)`
- Regressive: `a ∨ b` — ASCII: `CliffordAlgebras.regressiveprod(a,b)`
- Fat dot: `a ⋅ b` — ASCII: `CliffordAlgebras.fatdotprod(a,b)`
- Scalar: `a ⋆ b` — ASCII: `CliffordAlgebras.scalarprod(a,b)`
- Left contraction: `a ⨼ b` — ASCII: `CliffordAlgebras.leftcontractionprod(a,b)`
- Right contraction: `a ⨽ b` — ASCII: `CliffordAlgebras.rightcontractionprod(a,b)`
- Commutator: `a ×₋ b` — ASCII: `CliffordAlgebras.commutatorprod(a,b)`
- Anti-commutator: `a ×₊ b` — ASCII: `CliffordAlgebras.anticommutatorprod(a,b)`
- Sandwich (conjugation): `A ≀ X` — ASCII: `CliffordAlgebras.sandwichproduct(A,X)`

## Working with multivectors

```jldoctest cheatsheet
julia> cl3 = CliffordAlgebra(3);

julia> e1, e2, e3 = cl3.e1, cl3.e2, cl3.e3;

julia> mv = 1 + 2e1 + 3e2 + (e1 ∧ e2);  # 0-,1-,2- grades

julia> scalar(mv); grade(mv,2); even(mv); odd(mv); ~mv; dual(mv); true
true
```

## Rotors and motors (quick)

```jldoctest cheatsheet
julia> cl3 = CliffordAlgebra(3);

julia> B = (π/4) * (cl3.e1 ∧ cl3.e2); R = exp(B);

julia> v = cl3.e1; v_rot = R ≀ v; true
true
```

PGA motor (rotation+translation):

```jldoctest cheatsheet
julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3, e0 = basevector(pga,1), basevector(pga,2), basevector(pga,3), basevector(pga,:e0);

julia> B = (π/12)*(e1 ∧ e2); T = 0.05*(e0 ∧ e1); M = exp(B+T);

julia> v = e1 + 2e2; v2 = M ≀ v; true
true
```

## Outermorphism (linear transform)

```jldoctest cheatsheet
julia> cl2 = CliffordAlgebra(2); v = cl2.e1 + cl2.e2;

julia> θ = π/6; R = [cos(θ) -sin(θ); sin(θ) cos(θ)];

julia> v2 = outermorphism(R, v); true
true
```

## Tables and signatures

```jldoctest cheatsheet
julia> cl2 = CliffordAlgebra(2);

julia> io = IOBuffer(); signaturetable(io, cl2); cayleytable(io, cl2); true
true
```

## Quick facts

- Algebra dimension for signature (p,q,r) with n=p+q+r is 2^n.
- Number of grade-k elements is C(n,k).
- Pseudoscalar I has grade n; I^2 = character(algebra).
- Null basis vectors (r>0) square to 0 and represent infinity/conformal components depending on algebra.
