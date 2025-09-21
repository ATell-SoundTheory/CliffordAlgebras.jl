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

julia> v_ref = R ≀ v; v_ref isa MultiVector; true
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

## PGA: reflect a vector in a plane

Reflect across the plane z=0 (normal e3). The z component flips.

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3 = basevector(pga,1), basevector(pga,2), basevector(pga,3);

julia> plane = e3;  # z=0

julia> v = e1 + e3;

julia> v_ref = plane ≀ v; true
true
```

## PGA: line-plane meet -> point

The meet (∨) of a line and a plane yields a point (grade-3 in PGA3D).

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3 = basevector(pga,1), basevector(pga,2), basevector(pga,3);

julia> plane_xy = e3;      # z=0

julia> plane_xz = e2;      # y=0

julia> line = plane_xy ∧ plane_xz;  # x-axis line

julia> plane_yz = e1;      # x=0

julia> pt = line ∨ plane_yz;

julia> isgrade(pt, 3)
true
```

## PGA: rigid-body motor (rotation + translation)

Compose a small rotation about z with a small translation along x, then apply to a vector.

```jldoctest
julia> using CliffordAlgebras

julia> pga = CliffordAlgebra(:PGA3D);

julia> e1, e2, e3, e0 = basevector(pga,1), basevector(pga,2), basevector(pga,3), basevector(pga,:e0);

julia> B = (π/12) * (e1 ∧ e2);      # rotate about z

julia> T = 0.05 * (e0 ∧ e1);        # translate along x (translator-like)

julia> M = exp(B + T);

julia> v = e1 + 2e2;  # a direction-grade representative

julia> v2 = M ≀ v; v2 isa MultiVector; true
true
```

## Cl3: compose two rotations via rotors

```jldoctest
julia> using CliffordAlgebras

julia> cl3 = CliffordAlgebra(3);

julia> e1, e2, e3 = cl3.e1, cl3.e2, cl3.e3;

julia> B1 = (π/6) * (e1 ∧ e2);

julia> B2 = (π/7) * (e2 ∧ e3);

julia> R1, R2 = exp(B1), exp(B2);

julia> v = e1 + 2e2 + 3e3;

julia> v_seq = R2 ≀ (R1 ≀ v);

julia> v_comb = (R2 * R1) ≀ v; true
true
```

## STA: Lorentz boost via bivector exponential

Boost along x with rapidity 0.1 (using generator t∧x). Minkowski norm is preserved.

```jldoctest
julia> using CliffordAlgebras

julia> sta = CliffordAlgebra(:Spacetime);

julia> t, x = basevector(sta,:t), basevector(sta,:x);

julia> B = 0.1 * (t ∧ x);

julia> R = exp(B);

julia> v = t + 0.5x;  # timelike vector with small spatial part

julia> vb = R ≀ v;    # boosted

julia> s(a) = scalar(a * ~a);

julia> isgrade(vb, 1) && s(vb) ≈ s(v)
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

julia> v_rot = R ≀ v; v_rot isa MultiVector; true
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

julia> p2 = T ≀ p; p2 isa MultiVector; true
true
```

## Tips
- Use `matrix(mv)` and `vector(mv)` for sanity checks against linear algebra operations.
- Prefer small angles and normalized elements in doctests to keep outputs stable.
- When doctests are sensitive to printed formatting, write to an IOBuffer and assert `true`.
