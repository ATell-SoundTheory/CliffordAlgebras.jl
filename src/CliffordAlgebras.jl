module CliffordAlgebras

export CliffordAlgebra
export order, signature, dimension, character
export basesymbol, show 

export MultiVector
export zero, iszero, one, isone
export ==, algebra, eltype
export baseindices, coefficients
export scalar, prune, extend, grade
export grin, dual, even, odd, conj
export basevector, pseudoscalar
export vector, matrix

export +,-,*,/,\
export ∧, ∨, ⋅, ⨼, ⨽, ⋆
export ×₊, ×₋, ≀
export inv, conj, adjoint
export polarize, norm
export Λᵏ, exp

include("utils.jl")
include("algebra.jl")
include("multivector.jl")
include("arithmetic.jl")

end
