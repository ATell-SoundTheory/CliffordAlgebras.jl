# Tutorial: Getting Started with CliffordAlgebras.jl

This tutorial will walk you through the basics of using CliffordAlgebras.jl for geometric algebra computations.

## Installation

```jldoctest
julia> using CliffordAlgebras
```

## Lesson 1: Creating Your First Algebra

Let's start with the simplest non-trivial Clifford algebra, Cl(2,0,0):

```jldoctest
julia> cl2 = CliffordAlgebra(2)
Cl(2, 0, 0)

julia> io = IOBuffer(); cayleytable(io, cl2); true
true
```

This creates an algebra with 2 basis vectors e₁ and e₂ that both square to +1.

## Lesson 2: Basis Elements

Every Clifford algebra has 2ⁿ basis elements for n generators:

```jldoctest
julia> scalar_unit, e1, e2, e12 = cl2.𝟏, cl2.e1, cl2.e2, cl2.e1e2;

julia> e1*e1
1

julia> e2*e2
1

julia> e1*e2 == e12
true

julia> e2*e1 == -e12
true
```

## Lesson 3: Building Multivectors

A general multivector combines elements of different grades:

```jldoctest
julia> mv = 2.0 + 3.0*e1 + 4.0*e2 + 5.0*e12;

julia> scalar(mv)
2.0

julia> grade(mv, 1) isa MultiVector
true

julia> grade(mv, 2) isa MultiVector
true
```

## Lesson 4: The Geometric Product

The geometric product is the fundamental operation:

```jldoctest
julia> a = 2*e1 + 3*e2; b = 4*e1 + 5*e2;

julia> result = a * b; (scalar(result), grade(result, 2) isa MultiVector)
(23.0, true)
```

## Lesson 5: Specialized Products

### Exterior Product (Wedge Product)
Creates higher-grade elements:

```jldoctest
julia> area_element = e1 ∧ e2; area_element == e12
true

julia> e1 ∧ e1
0
```

### Interior Products
Various ways to contract multivectors:

```jldoctest
julia> fat_dot = a ⋅ b; fat_dot isa MultiVector
true

julia> scalar_prod = a ⋆ b; scalar(scalar_prod) isa Real
true
```

## Lesson 6: Rotations in 2D

One of the most important applications is representing rotations:

```jldoctest
julia> v = e1; angle = π/4; B = angle * e12;

julia> rotor = exp(B);

julia> rotated_v = rotor ≀ v;

julia> expected = (e1 + e2) / sqrt(2);

julia> rotated_v ≈ expected
true
```

## Lesson 7: Working with 3D Space

Let's move to three dimensions:

```julia
cl3 = CliffordAlgebra(3)

# Create some 3D vectors
x_axis = cl3.e1
y_axis = cl3.e2  
z_axis = cl3.e3

# Create a general vector
v3d = 1*x_axis + 2*y_axis + 3*z_axis

# The pseudoscalar in 3D
I3 = cl3.e1e2e3
println("3D pseudoscalar: ", I3)
println("I³² = ", I3 * I3)  # Should be -1

# Duality operation
dual_v = dual(v3d)
println("Dual of vector: ", dual_v)
```

## Lesson 8: Spacetime Algebra

Spacetime algebra uses signature (1,3,0):

```julia
sta = CliffordAlgebra(:Spacetime)

# Basis vectors
t = sta.t  # Time (squares to +1)
x = sta.x  # Space (squares to -1) 
y = sta.y
z = sta.z

# Create a spacetime event
event = 5*t + 3*x + 4*y + 0*z

# Compute the spacetime interval
interval = scalar(event * reverse(event))
println("Spacetime interval: ", interval)

# Check if timelike, spacelike, or lightlike
if interval > 0
    println("Timelike separation")
elseif interval < 0
    println("Spacelike separation")  
else
    println("Lightlike separation")
end
```

## Lesson 9: Advanced Operations

### Inverse and Division

```julia
cl2 = CliffordAlgebra(2)
mv = 2 + 3*cl2.e1

# Compute inverse
inv_mv = inv(mv)
println("Inverse: ", inv_mv)

# Verify: mv * inv(mv) should be 1
println("mv * inv(mv) = ", mv * inv_mv)

# Division
a = 1 + cl2.e1
b = 2 + cl2.e2
quotient = a / b
println("a / b = ", quotient)
```

### Exponential Function

```julia
# Exponential of bivectors gives rotors
B = π/6 * cl2.e1e2
rotor = exp(B)

# Exponential of vectors in hyperbolic case
hyp = CliffordAlgebra(1, 1)  # One positive, one negative
boost_generator = hyp.e1e2
boost = exp(0.5 * boost_generator)  # Hyperbolic rotation
```

## Next Steps

- Explore conformal geometric algebra for geometric modeling
- Learn about projective geometric algebra for computer graphics
- Study spacetime algebra for relativistic physics
- Investigate the package's performance features for large computations

## Common Pitfalls

1. **Forgetting anticommutativity**: Remember that `e1 * e2 = -e2 * e1`
2. **Mixing algebras**: You can't directly combine multivectors from different algebras
3. **Type stability**: Use `prune()` carefully as it's not type-stable
4. **Grade confusion**: Remember that the geometric product of two vectors has both grade 0 and grade 2 parts

## Performance Tips

1. Use type-stable operations when possible
2. Leverage the sparse representation for efficiency
3. Pre-compute rotors for repeated transformations
4. Use `@inferred` to check type stability in performance-critical code