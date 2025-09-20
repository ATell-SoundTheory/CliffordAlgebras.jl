[![Tests](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl/actions/workflows/Runtests.yml)
[![Documentation](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl/actions/workflows/Documenter.yml)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://atell-soundtheory.github.io/CliffordAlgebras.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://atell-soundtheory.github.io/CliffordAlgebras.jl/dev)


# CliffordAlgebras.jl

A fast, lightweight and easy-to-use Julia package for Clifford Algebras and Geometric Algebras.
Contributions welcome—see `CONTRIBUTING.md`. See `CHANGELOG.md` for notable changes. Many examples in the docs are doctested.
If you're contributing, see the Developer Guide in the docs: https://atell-soundtheory.github.io/CliffordAlgebras.jl/dev/developer/

CliffordAlgebras.jl uses compile time code generation, sparse multivector representation and special case identification to handle arbitrary geometric algebras efficiently. Lazy evaluation of expressions is not yet supported.

CliffordAlgebras provides a low level implementation that is common to all Clifford algebras. Functions that are specific to certain algebras can be added on top using the tools provided here.

Generate a Clifford algebra:

```
julia> cl2 = CliffordAlgebra(2)
Cl(2,0,0)

julia> cayleytable(stdout, cl2)
┌───────┬──────────────┬───────┐
│  +1   │  +e1    +e2  │ +e1e2 │
├───────┼──────────────┼───────┤
│  +e1  │  +1    +e1e2 │  +e2  │
│  +e2  │ -e1e2   +1   │  -e1  │
├───────┼──────────────┼───────┤
│ +e1e2 │  -e2    +e1  │  -1   │
└───────┴──────────────┴───────┘

julia> cl101 = CliffordAlgebra(1,0,1)
Cl(1,0,1)

julia> cayleytable(stdout, cl101)
┌───────┬──────────────┬───────┐
│  +1   │  +e1    +e2  │ +e1e2 │
├───────┼──────────────┼───────┤
│  +e1  │  +1    +e1e2 │  +e2  │
│  +e2  │ -e1e2    0   │   0   │
├───────┼──────────────┼───────┤
│ +e1e2 │  -e2     0   │   0   │
└───────┴──────────────┴───────┘

julia> sta = CliffordAlgebra(:Spacetime)
Cl(1,3,0)

julia> signaturetable(stdout, sta)
┌───┬────┐
│ t │ +1 │
│ x │ -1 │
│ y │ -1 │
│ z │ -1 │
└───┴────┘

julia> cayleytable(stdout, sta)
┌───────┬────────────────────────────┬──────────────────────────────────────────┬────────────────────────────┬───────┐
│  +1   │  +t     +x     +y     +z   │  +tx    +ty    +xy    +tz    +zx    +yz  │ +tyx   +txz   +tzy   +xyz  │ +txyz │
├───────┼────────────────────────────┼──────────────────────────────────────────┼────────────────────────────┼───────┤
│  +t   │  +1     +tx    +ty    +tz  │  +x     +y    -tyx    +z    -txz   -tzy  │  -xy    -zx    -yz   +txyz │ +xyz  │
│  +x   │  -tx    -1     +xy    -zx  │  +t    +tyx    -y    -txz    +z    +xyz  │  -ty    +tz   +txyz   -yz  │ -tzy  │
│  +y   │  -ty    -xy    -1     +yz  │ -tyx    +t     +x    +tzy   +xyz    -z   │  +tx   +txyz   -tz    -zx  │ -txz  │
│  +z   │  -tz    +zx    -yz    -1   │ +txz   -tzy   +xyz    +t     -x     +y   │ +txyz   -tx    +ty    -xy  │ -tyx  │
├───────┼────────────────────────────┼──────────────────────────────────────────┼────────────────────────────┼───────┤
│  +tx  │  -x     -t    -tyx   +txz  │  +1     -xy    -ty    +zx    +tz   +txyz │  -y     +z    +xyz   +tzy  │  +yz  │
│  +ty  │  -y    +tyx    -t    -tzy  │  +xy    +1     +tx    -yz   +txyz   -tz  │  +x    +xyz    -z    +txz  │  +zx  │
│  +xy  │ -tyx    +y     -x    +xyz  │  +ty    -tx    -1    +txyz   -yz    +zx  │  +t    -tzy   +txz    -z   │  -tz  │
│  +tz  │  -z    -txz   +tzy    -t   │  -zx    +yz   +txyz   +1     -tx    +ty  │ +xyz    -x     +y    +tyx  │  +xy  │
│  +zx  │ -txz    -z    +xyz    +x   │  -tz   +txyz   +yz    +tx    -1     -xy  │ +tzy    +t    -tyx    -y   │  -ty  │
│  +yz  │ -tzy   +xyz    +z     -y   │ +txyz   +tz    -zx    -ty    +xy    -1   │ -txz   +tyx    +t     -x   │  -tx  │
├───────┼────────────────────────────┼──────────────────────────────────────────┼────────────────────────────┼───────┤
│ +tyx  │  -xy    -ty    +tx   -txyz │  -y     +x     +t    -xyz   -tzy   +txz  │  -1     -yz    +zx    +tz  │  +z   │
│ +txz  │  -zx    +tz   -txyz   -tx  │  +z    -xyz   +tzy    -x     +t    -tyx  │  +yz    -1     -xy    +ty  │  +y   │
│ +tzy  │  -yz   -txyz   -tz    +ty  │ -xyz    -z    -txz    +y    +tyx    +t   │  -zx    +xy    -1     +tx  │  +x   │
│ +xyz  │ -txyz   -yz    -zx    -xy  │ -tzy   -txz    -z    -tyx    -y     -x   │  -tz    -ty    -tx    +1   │  -t   │
├───────┼────────────────────────────┼──────────────────────────────────────────┼────────────────────────────┼───────┤
│ +txyz │ -xyz   +tzy   +txz   +tyx  │  +yz    +zx    -tz    +xy    -ty    -tx  │  -z     -y     -x     +t   │  -1   │
└───────┴────────────────────────────┴──────────────────────────────────────────┴────────────────────────────┴───────┘
```

Get the basis vector names from the generated algebra and use the base vectors to construct multivectors:

```
julia> propertynames(cl2)
(:𝟏, :e1, :e2, :e1e2)

julia> mv1 = 2.0 * cl2.e1+ cl2.e1e2 + 1
+1.0+2.0×e1+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv2 = mv1 - cl2.e1e2
+1.0+2.0×e1 ∈ Cl(2, 0, 0)
```

Use the geometric product and derived products with the multivectors:
```
julia> mv1 * mv2
+5.0+4.0×e1-2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1 ∧ mv2
+1.0+4.0×e1+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1 ∨ mv2
+1.0+2.0×e1 ∈ Cl(2, 0, 0)

julia> mv1 ⋅ mv2
+5.0+4.0×e1-2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1 ⋆ mv2
+5.0 ∈ Cl(2, 0, 0)

julia> mv1 ⨼ mv2
+5.0+2.0×e1 ∈ Cl(2, 0, 0)

julia> mv1 ⨽ mv2
+5.0+2.0×e1-2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1 ×₋ mv2
-2.0×e2 ∈ Cl(2, 0, 0)

julia> mv1 ×₊ mv2
+5.0+4.0×e1+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1 ≀ mv2
+14.0+12.0×e1-8.0×e2 ∈ Cl(2, 0, 0)
```

Calculate the inverse of a multivector and divide multivectors:

``` 
julia> inv(mv2)
-0.33333333333333326+0.6666666666666665×e1 ∈ Cl(2, 0, 0)

julia> mv1 / mv2
+0.9999999999999998-0.6666666666666665×e2-0.33333333333333326×e1e2 ∈ Cl(2, 0, 0)

julia> mv2 \ mv1
+0.9999999999999998+0.6666666666666665×e2-0.33333333333333326×e1e2 ∈ Cl(2, 0, 0)
``` 

Calculate duals and involutions:

```
julia> dual(mv1)
+1.0+2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> grin(mv1)
+1.0-2.0×e1+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> conj(mv1)
+1.0-2.0×e1-1.0×e1e2 ∈ Cl(2, 0, 0)

julia> reverse(mv1)
+1.0+2.0×e1-1.0×e1e2 ∈ Cl(2, 0, 0)

julia> ~mv1
+1.0+2.0×e1-1.0×e1e2 ∈ Cl(2, 0, 0)

julia> polarize(mv1)
-1.0+2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> mv1'
-1.0+2.0×e2+1.0×e1e2 ∈ Cl(2, 0, 0)
```

Extract subspaces:

```
julia> scalar(mv1)
1.0

julia> even(mv1)
+1.0+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> odd(mv1)
+2.0×e1 ∈ Cl(2, 0, 0)

julia> grade(mv1,2)
+1.0×e1e2 ∈ Cl(2, 0, 0)

julia> Λᵏ(mv1,2)
+1.0×e1e2 ∈ Cl(2, 0, 0)
```

Calculate the norm of a multivector, evaluate the exponential function and use the sandwich product:

```
julia> R = cl2.e1e2 * π/8
+0.39269908169872414×e1e2 ∈ Cl(2, 0, 0)

julia> norm(R)
0.39269908169872414

julia> exp(R)
+0.9238795325112867+0.38268343236508984×e1e2 ∈ Cl(2, 0, 0)

julia> exp(R) ≀ cl2.e1
+0.7071067811865475×e1-0.7071067811865477×e2 ∈ Cl(2, 0, 0)
```

Map multivectors to their matrix algebra counterparts:

```
julia> vector(mv1)
4-element Vector{Float64}:
 1.0
 2.0
 0.0
 1.0

julia> matrix(mv1)
4×4 Matrix{Float64}:
 1.0   2.0  0.0  -1.0
 2.0   1.0  1.0   0.0
 0.0  -1.0  1.0   2.0
 1.0   0.0  2.0   1.0  

julia> matrix(mv1)*matrix(mv2) == matrix(mv1*mv2)
true

julia> matrix(mv2)*vector(mv1) == vector(mv2*mv1)
true
```

Apply outermorphisms:
```
julia> M = [0 1 ; 1 0]
2×2 Matrix{Int64}:
 0  1
 1  0

julia> outermorphism(M, mv1)
+1.0+2.0×e2-1.0×e1e2 ∈ Cl(2, 0, 0)
```
