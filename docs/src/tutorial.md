# Tutorial: Getting Started with CliffordAlgebras.jl

This tutorial will walk you through the basics of using CliffordAlgebras.jl for geometric algebra computations.

## Installation

```jldoctest
using Pkg
Pkg.add("CliffordAlgebras")
using CliffordAlgebras
```

## Lesson 1: Creating Your First Algebra

Let's start with the simplest non-trivial Clifford algebra, Cl(2,0,0):

```jldoctest
# Create a 2D Euclidean geometric algebra
cl2 = CliffordAlgebra(2)

# Inspect the algebra
println(cl2)  # Prints: Cl(2,0,0)

# View the multiplication table (smoke test)
cayleytable(IOBuffer(), cl2)
true
```

This creates an algebra with 2 basis vectors e₁ and e₂ that both square to +1.

## Lesson 2: Basis Elements

Every Clifford algebra has 2ⁿ basis elements for n generators:

```jldoctest
# Access basis elements
scalar_unit = cl2.𝟏      # Scalar unit
e1 = cl2.e1             # First basis vector  
e2 = cl2.e2             # Second basis vector
e12 = cl2.e1e2          # Bivector (area element)

# Check their properties
println("e1² = ", e1 * e1)     # Should be +1
println("e2² = ", e2 * e2)     # Should be +1  
println("e1e2 = ", e1 * e2)    # Should be e1e2
println("e2e1 = ", e2 * e1)    # Should be -e1e2 (anticommutative)
```

## Lesson 3: Building Multivectors

A general multivector combines elements of different grades:

```jldoctest
# Create a general multivector
mv = 2.0 + 3.0*e1 + 4.0*e2 + 5.0*e12

println("Multivector: ", mv)

# Extract components by grade
println("Scalar part: ", scalar(mv))
println("Vector part: ", grade(mv, 1))
println("Bivector part: ", grade(mv, 2))
```

## Lesson 4: The Geometric Product

The geometric product is the fundamental operation:

```jldoctest
# Two vectors
a = 2*e1 + 3*e2
b = 4*e1 + 5*e2

# Geometric product
result = a * b
println("a * b = ", result)

# The result has both scalar and bivector parts
println("Scalar part: ", scalar(result))        # Dot product
println("Bivector part: ", grade(result, 2))    # Wedge product
```

## Lesson 5: Specialized Products

### Exterior Product (Wedge Product)
Creates higher-grade elements:

```jldoctest
# Exterior product of vectors creates bivector
area_element = e1 ∧ e2
println("e1 ∧ e2 = ", area_element)

# Vectors wedge with themselves give zero
println("e1 ∧ e1 = ", e1 ∧ e1)
```

### Interior Products
Various ways to contract multivectors:

```jldoctest
# Fat dot product (symmetric part of geometric product)
fat_dot = a ⋅ b
println("a ⋅ b = ", fat_dot)

# Scalar product (fully contracted)
scalar_prod = a ⋆ b  
println("a ⋆ b = ", scalar_prod)
```

## Lesson 6: Rotations in 2D

One of the most important applications is representing rotations:

```jldoctest
# Create a vector to rotate
v = e1

# Create a rotation bivector (angle in bivector form)
angle = π/4  # 45 degrees
B = angle * e12

# Create the rotor using exponential
rotor = exp(B)
println("Rotor: ", rotor)

# Apply rotation using sandwich product
rotated_v = rotor ≀ v
println("Rotated vector: ", rotated_v)

# Verify the rotation (should be (e1 + e2)/√2)
expected = (e1 + e2) / sqrt(2)
println("Expected: ", expected)
println("Close? ", rotated_v ≈ expected)
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