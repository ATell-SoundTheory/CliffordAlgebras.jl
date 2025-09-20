# Performance Guide

This guide provides tips and best practices for achieving optimal performance with CliffordAlgebras.jl.

## Understanding the Design

CliffordAlgebras.jl is designed for high performance through several key features:

1. **Sparse Representation**: Only non-zero coefficients are stored
2. **Generated Functions**: Compile-time specialization for each algebra
3. **Type Stability**: Most operations have inferrable return types
4. **Zero-Cost Abstractions**: Overhead is eliminated at compile time

## Type Stability

Type stability is crucial for Julia performance. Most CliffordAlgebras.jl operations are type-stable:

```julia
using CliffordAlgebras

cl3 = CliffordAlgebra(3)
mv1 = cl3.e1 + cl3.e2
mv2 = cl3.e2 + cl3.e3

# These operations are type-stable
@inferred mv1 * mv2          # Geometric product
@inferred mv1 ∧ mv2          # Exterior product  
@inferred ~mv1               # Reverse
@inferred exp(cl3.e1e2)      # Exponential
```

### Non-Type-Stable Operations

Some operations are intentionally not type-stable for flexibility:

```julia
# prune() changes the sparse structure
pruned = prune(mv1)  # Return type depends on which coefficients survive

# Use extend() for type stability if needed
extended = extend(mv1)  # Always returns full representation
```

## Memory Efficiency

### Sparse vs Dense Representation

CliffordAlgebras.jl automatically uses sparse representation:

```julia
cl8 = CliffordAlgebra(8)  # 2^8 = 256 basis elements

# This only stores 2 coefficients, not 256
sparse_mv = cl8.e1 + cl8.e8

# Convert to dense if needed (usually not recommended)
dense_mv = extend(sparse_mv)
```

### Memory Usage Guidelines

1. **Prefer sparse operations**: Most functions preserve sparsity
2. **Avoid unnecessary extend()**: Only use when type stability is critical
3. **Reuse multivectors**: Modify in-place when possible

## Compile-Time Optimization

### Generated Functions

The package uses `@generated` functions for optimal performance:

```julia
# This generates specialized code for each algebra and multivector type
function (*)(a::MultiVector{CA}, b::MultiVector{CA}) where CA
    # Specialized implementation generated at compile time
end
```

### Specialization Tips

1. **Use concrete types**: Avoid `MultiVector` without type parameters
2. **Consistent algebra types**: Don't mix algebras in hot loops
3. **Stable multivector structures**: Reuse similar sparse patterns

## Benchmarking Examples

### Basic Operations

```julia
using BenchmarkTools, CliffordAlgebras

cl3 = CliffordAlgebra(3)
mv1 = 1.0 + cl3.e1 + cl3.e2 + cl3.e1e2
mv2 = 2.0 + cl3.e2 + cl3.e3 + cl3.e2e3

# Benchmark geometric product
@benchmark $mv1 * $mv2

# Benchmark exterior product  
@benchmark $mv1 ∧ $mv2

# Benchmark exponential
B = π/4 * cl3.e1e2
@benchmark exp($B)
```

### Large Algebras

```julia
# For larger algebras, sparsity becomes crucial
cl6 = CliffordAlgebra(6)  # 64 basis elements

# Sparse multivector (only 3 coefficients)
sparse_mv = cl6.e1 + cl6.e3 + cl6.e6

# Still efficient due to sparsity
@benchmark $sparse_mv * $sparse_mv
```

## Optimization Strategies

### 1. Precompute Common Operations

```julia
# Instead of recomputing rotors
angle = π/6
B = angle * cl3.e1e2
rotor = exp(B)  # Precompute this

# Use the precomputed rotor many times
for vector in vectors
    rotated = rotor ≀ vector
end
```

### 2. Use Appropriate Signatures

Choose the minimal signature for your problem:

```julia
# For 2D rotations, Cl(2,0,0) is more efficient than Cl(3,0,0)
cl2 = CliffordAlgebra(2)  # 4 basis elements vs 8

# For spacetime, use the exact signature
sta = CliffordAlgebra(1, 3, 0)  # Not Cl(4,0,0)
```

### 3. Minimize Allocations

```julia
# Good: Reuse multivectors when possible
function rotate_many_vectors!(results, rotor, vectors)
    for (i, v) in enumerate(vectors)
        results[i] = rotor ≀ v
    end
end

# Avoid: Creating new algebras in hot loops
function bad_example()
    for i in 1:1000
        cl = CliffordAlgebra(3)  # Don't do this!
        # ... operations
    end
end
```

### 4. Leverage Type Annotations

```julia
# Help the compiler with type annotations
function efficient_computation(mv::MultiVector{CA,Float64}) where CA
    result = mv * mv
    return scalar(result)
end
```

## Performance Pitfalls

### 1. Type Instability

```julia
# Bad: Type-unstable function
function unstable_norm(mv)
    if some_condition
        return norm(mv)  # Returns Float64
    else
        return mv        # Returns MultiVector
    end
end

# Good: Type-stable alternatives
function stable_norm(mv)
    return norm(mv)  # Always returns Float64
end
```

### 2. Unnecessary Conversions

```julia
# Bad: Converting between algebras
cl2_mv = cl2.e1
cl3_mv = MultiVector(cl3, (1.0, 0.0, 0.0))  # Expensive conversion

# Good: Work within one algebra
mv1 = cl3.e1
mv2 = cl3.e2
result = mv1 * mv2
```

### 3. Overuse of `prune()`

```julia
# Bad: Excessive pruning
result = mv1 * mv2
result = prune(result)  # Type-unstable and often unnecessary
result = result + mv3
result = prune(result)  # Again!

# Good: Prune only when needed
result = mv1 * mv2 + mv3
# Only prune if you know there are many small coefficients
if need_cleanup
    result = prune(result)
end
```

## Profiling Tools

### Memory Allocation

```julia
using Profile

function profile_example()
    cl4 = CliffordAlgebra(4)
    mv = cl4.e1 + cl4.e2 + cl4.e3 + cl4.e4
    
    for i in 1:1000
        result = mv * mv
    end
end

# Profile memory allocations
@profile profile_example()
Profile.print()
```

### Type Inference

```julia
using Cthulhu

# Inspect generated code
cl3 = CliffordAlgebra(3)
mv1 = cl3.e1
mv2 = cl3.e2

# Descend into the multiplication
@descend mv1 * mv2
```

## Performance Summary

**Fast Operations**:
- Geometric product between similar multivectors
- Exponential of bivectors  
- Grade extraction
- Reverse and other involutions

**Moderate Operations**:
- Operations between very different sparse structures
- Converting between representations
- Complex expressions with many terms

**Slow Operations**:
- Excessive use of `prune()`
- Mixing different algebras
- Type-unstable code patterns

**Best Practices**:
1. Stick to one algebra type per computation
2. Leverage sparsity naturally
3. Precompute rotors and other expensive operations
4. Use `@inferred` to check type stability
5. Profile before optimizing
6. Choose minimal sufficient algebra signatures