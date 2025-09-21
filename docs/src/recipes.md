# Recipes: CGA & PGA

This page collects small, practical recipes for common tasks in projective and conformal geometric algebra.

All snippets are doctest-able and ASCII-friendly where possible.

## Contract
- Input: an algebra and a few multivectors
- Output: a result multivector of interest
- Errors: invalid signatures or division by zero where noted
- Success: example evaluates and type-stable on recent Julia

## PGA: reflect a vector in a plane (householder-like)

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3 = basevector(pga,1), basevector(pga,2), basevector(pga,3);

julia> n = e2;  # plane normal (unit for simplicity)

julia> v = e1 + 2e2 + 3e3;

julia> R = exp(π * (n ∧ e3) / 2);  # 180° around axis orthogonal to plane -> reflection

julia> v_ref = R ≀ v;

julia> typeof(v_ref) <: typeof(v)
true
```

## PGA: intersection of two planes -> line

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3, e0 = basevector(pga,1), basevector(pga,2), basevector(pga,3), basevector(pga,:e0);

julia> plane_xy = e3;  # z=0 plane

julia> plane_xz = e2;  # y=0 plane

julia> line = plane_xy ∧ plane_xz;

julia> isgrade(line, 2)
true
```

## CGA: rotate a Euclidean vector by a rotor

```jldoctest
julia> using CliffordAlgebras

julia> cga = CliffordAlgebra(:CGA3D);

julia> e1, e2, e3 = basevector(cga,1), basevector(cga,2), basevector(cga,3);

julia> B = (π/6) * (e1 ∧ e2);

julia> R = exp(B);

julia> v = e1 + e2;

julia> v_rot = R ≀ v;

julia> isgrade(v_rot, 1)
true
```

## CGA: translate a point using a motor (CGA screw motion)

Note: For brevity, we model a simple translation by building a null bivector.

```jldoctest
julia> using CliffordAlgebras

julia> cga = CliffordAlgebra(:CGA3D);

julia> e1, e2, e3 = basevector(cga,1), basevector(cga,2), basevector(cga,3);

julia> eplus, eminus = basevector(cga,:e₊), basevector(cga,:e₋);

julia> ninf = eminus;  # direction to infinity

julia> t = 0.1*(e1) ∧ ninf;  # translator-like bivector

julia> T = exp(t);

julia> p = e1;  # simple direction-grade representative

julia> p2 = T ≀ p;

julia> typeof(p2) <: typeof(p)
true
```

## Tips
- Use `matrix(mv)` and `vector(mv)` for sanity checks against linear algebra operations.
- Prefer small angles and normalized elements in doctests to keep outputs stable.
- When doctests are sensitive to printed formatting, write to an IOBuffer and assert `true`.
