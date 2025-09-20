module CliffordAlgebras

export CliffordAlgebra
export order, signature, dimension, character
export basesymbol, show, show_multivector
export cayleytable, signaturetable

export MultiVector
export zero, iszero, one, isone
export ==, algebra, eltype, isapprox
export baseindices, coefficients, coefficient
export scalar, prune, extend
export grade, isgrade, maxgrade, mingrade
export reverse, ~
export conj, grin, dual, even, odd
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
