# CliffordAlgebras.jl Documentation

A fast, lightweight and easy-to-use Julia package for Clifford Algebras and Geometric Algebras.

## Overview

CliffordAlgebras.jl provides a comprehensive implementation of Clifford algebras (also known as geometric algebras) in Julia. The package uses compile-time code generation, sparse multivector representation, and special case identification to handle arbitrary geometric algebras efficiently.

## Quick Start

Install and load the package:

```julia
import Pkg
Pkg.add(url="https://github.com/ATell-SoundTheory/CliffordAlgebras.jl")
using CliffordAlgebras
```

```jldoctest
julia> using CliffordAlgebras

julia> cl2 = CliffordAlgebra(2)
Cl(2,0,0)

julia> e1, e2, e12 = cl2.e1, cl2.e2, cl2.e1e2;

julia> mv = 1 + 2*e1 + 3*e2 + 4*e12;

julia> mv * e1 isa MultiVector
true

julia> e1 ∧ e2 == cl2.e1e2
true
```

## Mathematical Background

### Clifford Algebras

A Clifford algebra Cl(p,q,r) is an associative algebra over the real numbers with a signature (p,q,r), where:
- p = number of basis vectors that square to +1
- q = number of basis vectors that square to -1  
- r = number of basis vectors that square to 0

The fundamental relation is the Clifford relation:
```
eᵢeⱼ + eⱼeᵢ = 2ηᵢⱼ
```
where ηᵢⱼ is the metric tensor.

### Multivectors

A general multivector in Cl(p,q,r) can be written as:
```
M = α + aᵢeᵢ + aᵢⱼeᵢeⱼ + ... + a₁₂...ₙe₁e₂...eₙ
```

The package represents multivectors using a sparse representation, storing only non-zero coefficients.

## Creating Algebras

### Basic Signatures

```julia
# Euclidean algebras
cl2 = CliffordAlgebra(2)        # Cl(2,0,0)
cl3 = CliffordAlgebra(3)        # Cl(3,0,0)

# Mixed signatures
spa = CliffordAlgebra(1,3)      # Cl(1,3,0) - Spacetime algebra
pga = CliffordAlgebra(3,0,1)    # Cl(3,0,1) - Projective geometric algebra

# Full signature specification
custom = CliffordAlgebra(2,1,1) # Cl(2,1,1)
```

### Predefined Algebras

The package includes many commonly used algebras:

```julia
# Complex numbers
ℂ = CliffordAlgebra(:Complex)

# Quaternions  
ℍ = CliffordAlgebra(:Quaternions)

# Spacetime algebra
sta = CliffordAlgebra(:Spacetime)

# Projective geometric algebras
pga2d = CliffordAlgebra(:PGA2D)
pga3d = CliffordAlgebra(:PGA3D)

# Conformal geometric algebras
cga2d = CliffordAlgebra(:CGA2D)
cga3d = CliffordAlgebra(:CGA3D)
```

Tip: See “Typing Unicode operators” in the README for how to enter symbols like ∧, ⋅, ⨼, and ≀ in the REPL and editors.

## Working with Multivectors

### Creating Multivectors

```julia
cl3 = CliffordAlgebra(3)

# From scalar
scalar_mv = MultiVector(cl3, 5.0)

# From basis vectors
mv1 = cl3.e1 + 2*cl3.e2 + 3*cl3.e1e2

# From vector components
vector_mv = MultiVector(cl3, (1.0, 2.0, 3.0))
```

### Operations

#### Geometric Product (*)
The fundamental operation in Clifford algebra:
```julia
result = mv1 * mv2
```

#### Exterior Product (∧)
Anti-commutative product that creates higher-grade elements:
```julia
bivector = cl3.e1 ∧ cl3.e2  # Creates e1e2
```

If you prefer ASCII, you can use function names instead of Unicode operators, e.g. `exteriorprod(mv1, mv2)` for `mv1 ∧ mv2`.

#### Interior Products
Various contraction operations:
```julia
left_contraction = mv1 ⨼ mv2   # Left contraction
right_contraction = mv1 ⨽ mv2  # Right contraction
fat_dot = mv1 ⋅ mv2            # Fat dot product
scalar_product = mv1 ⋆ mv2     # Scalar product
```

#### Involutions
```julia
reverse_mv = ~mv1              # Reverse (grade involution)
conjugate_mv = conj(mv1)       # Clifford conjugation
grade_inv = grin(mv1)          # Grade involution
dual_mv = dual(mv1)            # Hodge dual
```

### Grade Operations

```julia
# Extract specific grades
scalar_part = scalar(mv)
grade_2 = grade(mv, 2)
even_part = even(mv)
odd_part = odd(mv)

# Grade queries
is_vector = isgrade(mv, 1)
max_grade_mv, max_g = maxgrade(mv)
```

## Advanced Features

### Exponential and Trigonometric Functions

```julia
# Exponential of bivector (rotation)
B = π/4 * cl3.e1e2
rotor = exp(B)

# Apply rotation using sandwich product
rotated_vector = rotor ≀ cl3.e3
```

### Matrix Representation

```julia
# Convert to matrix form
matrix_form = matrix(mv)
vector_form = vector(mv)

# Matrix multiplication is equivalent to geometric product
@assert matrix(mv1) * matrix(mv2) == matrix(mv1 * mv2)
```

### Outermorphisms

Apply linear transformations to the underlying vector space:

```julia
# 2x2 rotation matrix
θ = π/4
R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

# Apply to multivector
transformed_mv = outermorphism(R, mv)
```

## Performance Considerations

### Sparse Representation
The package automatically uses sparse representation for efficiency:
```julia
# Only stores non-zero coefficients
sparse_mv = cl3.e1 + cl3.e1e2e3  # Only 2 coefficients stored
```

### Type Stability
Most operations are type-stable for optimal performance:
```julia
# These operations are typically type-stable
@inferred mv1 * mv2
@inferred exp(bivector)
@inferred ~mv1
```

### Generated Functions
The package uses `@generated` functions for compile-time optimization, ensuring efficient code for each algebra and multivector type.

## Examples

### 2D Rotations
```julia
cl2 = CliffordAlgebra(2)

# Create a vector
v = cl2.e1

# Create rotation bivector (π/4 rotation)
B = π/4 * cl2.e1e2

# Create rotor
R = exp(B)

# Apply rotation
rotated_v = R ≀ v
```

### 3D Spacetime
```julia
sta = CliffordAlgebra(:Spacetime)  # Cl(1,3,0)

# Create spacetime vector
event = sta.t + 2*sta.x + 3*sta.y + 4*sta.z

# Compute spacetime interval
interval = scalar(event * ~event)
```

 
