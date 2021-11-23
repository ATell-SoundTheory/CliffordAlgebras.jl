module CliffordAlgebras

export CliffordAlgebra
export order, signature, dimension, character
export basesymbol, show

export MultiVector
export zero, iszero, one, isone
export ==, algebra, eltype
export baseindices, coefficients
export scalar, prune, extend
export grade, reverse, ~, conj
export grin, dual, even, odd
export basevector, pseudoscalar
export vector, matrix

export +, -, *, /, \
export ∧, ∨, ⋅, ⨼, ⨽, ⋆
export ×₊, ×₋, ≀
export inv, adjoint
export polarize, norm, norm_sqr
export Λᵏ, exp, outermorphism

include("utils.jl")
include("algebra.jl")
include("multivector.jl")
include("arithmetic.jl")

end
