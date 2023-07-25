# MultiVector arithmetic for CliffordAlgebras.jl

import Base.+, Base.-, Base.*, Base./, Base.\
import Base.broadcasted
import Base.inv, Base.adjoint, Base.exp
import LinearAlgebra.norm, LinearAlgebra.norm_sqr
import LinearAlgebra.SingularException
import StaticArrays.SVector, StaticArrays.SMatrix
import SparseArrays.sparse, SparseArrays.sparsevec


@generated function (+)(a::MultiVector{CA,Ta}, b::MultiVector{CA,Tb}) where {CA,Ta,Tb}
    d = dimension(CA)
    acc_a = [Int(0) for n = 1:d]
    acc_b = [Int(0) for n = 1:d]
    bi_a = baseindices(a)
    bi_b = baseindices(b)
    acc_a[collect(bi_a)] = 1:length(bi_a)
    acc_b[collect(bi_b)] = 1:length(bi_b)
    BI = Tuple(findall(map( (a,b) -> !iszero(a) || !iszero(b), acc_a, acc_b)))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(Algebra(a), zero($T)))
    end
    coeffs = (
        iszero(acc_a[n]) ? :(coefficients(b)[$(acc_b[n])]) :
        iszero(acc_b[n]) ? :(coefficients(a)[$(acc_a[n])]) : :(coefficients(a)[$(acc_a[n])] + coefficients(b)[$(acc_b[n])]) for
        n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(@inbounds MultiVector{$CA,$T,$BI,$K}($coeffsexpr))
end

(+)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA, a) + b
(+)(a::MultiVector{CA}, b::Real) where {CA} = a + MultiVector(CA, b)

(-)(a::MultiVector) = map_coefficients(-, a)
(-)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = a + (-b)
(-)(a::Real, b::MultiVector) = a + (-b)
(-)(a::MultiVector, b::Real) = a + (-b)


function generatefilteredprod(
    a::Type{<:MultiVector{CA,Ta}},
    b::Type{<:MultiVector{CA,Tb}},
    ::Val{Filter},
) where {CA,Ta,Tb,Filter}
    d = dimension(CA)
    acc = [NTuple{3,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)), (kb, base_b) in enumerate(baseindices(b))
        baseidx, coeff = baseproduct(CA, base_a, base_b)
        leftgrade = basegrade(CA, base_a)
        rightgrade = basegrade(CA, base_b)
        productgrade = basegrade(CA, baseidx)
        keep = if Filter == :geometric
            true
        elseif Filter == :exterior
            productgrade == leftgrade + rightgrade
        elseif Filter == :rightcontraction
            productgrade == leftgrade - rightgrade
        elseif Filter == :leftcontraction
            productgrade == rightgrade - leftgrade
        elseif Filter == :fatdot
            productgrade == abs(leftgrade - rightgrade)
        elseif Filter == :scalar
            iszero(productgrade)
        else
            error("unknown filter.")
        end
        if !iszero(coeff) && keep
            push!(acc[baseidx], (ka, kb, coeff))
        end
    end
    BI = Tuple(findall(!isempty, acc))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(Algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (:($(tup[3]) * coefficients(a)[$(tup[1])] * coefficients(b)[$(tup[2])]) for tup in acc[n])...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(@inbounds MultiVector{$CA,$T,$BI,$K}($coeffsexpr))
end

@generated function geometricprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:geometric))
end

@generated function exteriorprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:exterior))
end

@generated function leftcontractionprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:leftcontraction))
end

@generated function rightcontractionprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:rightcontraction))
end

@generated function fatdotprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:fatdot))
end

@generated function scalarprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    generatefilteredprod(a,b,Val(:scalar))
end


"""
    a * b

Calculates the geometric product of two MultiVectors a and b.
"""
(*)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = geometricprod(a, b)
(*)(a::Real, b::MultiVector) = mul_with_scalar(a,b)
(*)(a::MultiVector, b::Real) = mul_with_scalar(b,a)

function map_coefficients(f, mv::MultiVector{CA,T,BI,K}) where {CA,T,BI,K}
    new_coeffs = map(f, coefficients(mv))
    T_new = eltype(new_coeffs)
    MultiVector{CA,T_new, BI, K}(new_coeffs)
end

function mul_with_scalar(s::Real,mv::MultiVector)
    map_coefficients(x->s*x,mv)
end

"""
    mv1 .* mv2

Calculates the element wise product between the Mulivectors a and b.
"""
function broadcasted(::typeof(*), a::MultiVector{CA,Ta,BI}, b::MultiVector{CA,Tb,BI})::MultiVector where {CA,Ta,Tb,BI}
    return MultiVector(CA, BI, coefficients(a) .* coefficients(b))
end

function broadcasted(::typeof(*), a::MultiVector{CA,Ta,BIa}, b::MultiVector{CA,Tb,BIb})::MultiVector where {CA,Ta,Tb,BIa,BIb}
    BI = tuple(union(BIa, BIb)...)
    v1, v2 = coefficients(a, BI), coefficients(b, BI)
    return MultiVector(CA, BI, v1 .* v2)
end

"""
    a ∧ b

Calculates the wedge product between two MultiVectors a and b.
"""
(∧)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = exteriorprod(a, b)
(∧)(a::MultiVector{CA}, b::Real) where {CA} = mul_with_scalar(b,a)
(∧)(a::Real, b::MultiVector{CA}) where {CA} = mul_with_scalar(a,b)
(∧)(a::Real, b::Real) = a * b


"""
    a ⋅ b

Calculates the "fat dot" product between the MultiVectors a and b.
"""
(⋅)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = fatdotprod(a, b)
(⋅)(a::MultiVector{CA}, b::Real) where {CA} = mul_with_scalar(b,a)
(⋅)(a::Real, b::MultiVector{CA}) where {CA} = mul_with_scalar(a,b)
(⋅)(a::Real, b::Real) = a * b


"""
    a ⨼ b

Calculates the left contraction of the MultiVectors a and b.
"""
(⨼)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = leftcontractionprod(a, b)
(⨼)(a::MultiVector{CA}, b::Real) where {CA} = a ⨼ MultiVector(CA,b)
(⨼)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA,a) ⨼ b
(⨼)(a::Real, b::Real) = a * b


"""
    a ⨽ b

Calculates the right contraction of the MultiVectors a and b.
"""
(⨽)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = rightcontractionprod(a, b)
(⨽)(a::MultiVector{CA}, b::Real) where {CA} = a ⨽ MultiVector(CA,b)
(⨽)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA,a) ⨽ b
(⨽)(a::Real, b::Real) = a * b

"""
    a ⋆ b

Calculates the scalar product of the MultiVectors a and b.
"""
(⋆)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = scalarprod(a, b)
(⋆)(a::MultiVector{CA}, b::Real) where {CA} = a ⋆ MultiVector(CA,b)
(⋆)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA,a) ⋆ b
(⋆)(a::Real, b::Real) = a * b


"""
    a ∨ b

Calculates the regressive product of the MultiVectors a and b.
"""
(∨)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = dual(dual(a) ∧ dual(b))
(∨)(a::MultiVector{CA}, b::Real) where {CA} = a ∨ MultiVector(CA,b)
(∨)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA,a) ∨ b
(∨)(a::Real, b::Real) = zero(promote_type(typeof(a),typeof(b)))


function generatesymprod(
    a::Type{<:MultiVector{CA,Ta}},
    b::Type{<:MultiVector{CA,Tb}},
    ::Val{P}
) where {CA,Ta,Tb,P}
    d = dimension(CA)
    acc = [NTuple{3,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)), (kb, base_b) in enumerate(baseindices(b))
        baseidx1, coeff1 = baseproduct(CA, base_a, base_b)
        baseidx2, coeff2 = baseproduct(CA, base_b, base_a)
        if baseidx1 == baseidx2
            coeff = coeff1 + P * coeff2
            if !iszero(coeff)
                push!(acc[baseidx1], (ka, kb, coeff))
            end
        else
            if !iszero(coeff1)
                push!(acc[baseidx1], (ka, kb, coeff1))
            end
            if !iszero(coeff2)
                push!(acc[baseidx2], (ka, kb, P * coeff2))
            end
        end
    end
    BI = Tuple(findall(!isempty, acc))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(Algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (:($(tup[3]) * coefficients(a)[$(tup[1])] * coefficients(b)[$(tup[2])]) for tup in acc[n])...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(@inbounds MultiVector(Algebra(a), $BI, $coeffsexpr))
end


@generated function commutatorprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb}
) where {CA,Ta,Tb}
    generatesymprod(a, b, Val(-1))
end

@generated function anticommutatorprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb}
) where {CA,Ta,Tb}
    generatesymprod(a, b, Val(+1))
end


"""
    a ×₋ b

Calculates the commutator ab-ba of two MultiVectors a and b.
"""
(×₋)(a::MultiVector{CA,Ta}, b::MultiVector{CA,Tb}) where {CA,Ta,Tb} =
    inv(promote_type(Ta, Tb)(2)) * commutatorprod(a, b)

(×₋)(::Ta, ::MultiVector{CA,Tb}) where {CA,Ta<:Real,Tb} = zero(promote_type(Ta, Tb))
(×₋)(a::MultiVector, b::Real) = b ×₋ a
(×₋)(::Ta, ::Tb) where {Ta<:Real,Tb<:Real} = zero(promote_type(Ta, Tb))

"""
    a ×₊ b

Calculates the anti-commutator ab+ba of two MultiVectors a and b.
"""
(×₊)(a::MultiVector{CA,Ta}, b::MultiVector{CA,Tb}) where {CA,Ta,Tb} =
    inv(promote_type(Ta, Tb)(2)) * anticommutatorprod(a, b)
(×₊)(a::Real, b::MultiVector) = a * b
(×₊)(a::MultiVector, b::Real) = a * b
(×₊)(a::Real, b::Real) = a * b


@generated function sandwichproduct(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    # a * b * reverse(a)
    d = dimension(CA)
    acc = [NTuple{4,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)),
        (kb, base_b) in enumerate(baseindices(b)),
        (ka_rev, base_a_rev) in enumerate(baseindices(a))

        leftbaseidx, leftcoeff = baseproduct(CA, base_a, base_b)
        baseidx, rightcoeff = baseproduct(CA, leftbaseidx, base_a_rev)
        coeff = leftcoeff * rightcoeff * basereverse(CA, base_a_rev)
        if !iszero(coeff)
            push!(acc[baseidx], (ka, kb, ka_rev, coeff))
        end
    end
    BI = Tuple(findall(!isempty, acc))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(Algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (
                :($(tup[4]) * coefficients(a)[$(tup[1])] * coefficients(b)[$(tup[2])] * coefficients(a)[$(tup[3])]) for
                tup in acc[n]
            )...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(@inbounds MultiVector{$CA,$T,$BI,$K}($coeffsexpr))
end

"""
    a ≀ b

Calculates the sandwich product `a*b*reverse(a)` for two MultiVectors a and b.
"""
(≀)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = sandwichproduct(a, b)
(≀)(a::Real, b::MultiVector{CA}) where {CA} = a^2 * b
(≀)(a::MultiVector{CA}, b::Real) where {CA} = b * a *reverse(a)
(≀)(a::Real, b::Real) = a^2 * b


"""
    polarize(mv::MultiVector)
    mv'

Calculates the polarization of the MultiVector, i.e. `mv * pseudoscalar`.
"""
polarize(mv::MultiVector{CA}) where {CA} = mv * pseudoscalar(CA)

adjoint(mv::MultiVector) = polarize(mv)


"""
    Λᵏ(::MultiVector, ::Val{k}) where k
    Λᵏ(::MultiVector, k::Integer) 

Projects the MultiVector onto k-vectors. Similar to grade(mv,k), but uses
@generated code and compile time optimizations.
"""
@generated function (Λᵏ)(mv::MultiVector{CA}, ::Val{k}) where {CA,k}
    indexbounds =
        1 .+ (0, cumsum(ntuple(i -> binomial(order(CA), i - 1), 1 + order(CA)))...)
    @assert 0 < k + 1 < length(indexbounds) "k out of bounds."
    bi_low = indexbounds[k+1]
    bi_high = indexbounds[k+2] - 1
    s = findall(i -> i >= bi_low && i <= bi_high, baseindices(mv))
    BI = baseindices(mv)[s]
    K = length(BI)
    T = eltype(mv)
    coeffs = (:(coefficients(mv)[$i]) for i in s)
    if length(coeffs) == 0
        return :(MultiVector($CA, zero($T)))
    end
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(@inbounds MultiVector{$CA,$T,$BI,$K}($coeffsexpr))
end

(Λᵏ)(mv::MultiVector, k::Integer) = Λᵏ(mv, Val(k))

"""
    norm(::MultiVector)

Calculates the MultiVector norm defined as `sqrt(scalar(mv*reverse(mv)))`.
"""
function norm(mv::MultiVector{CA}) where CA
    Nn = signature(CA)[2]
    if Nn != 0
        sqrt(complex(norm_sqr(mv)))
    else
        sqrt(norm_sqr(mv))
    end
end

function default_atol(mv1, mv2)
    T = promote_type(eltype(mv1), eltype(mv2))
    zero(T)
end
function default_rtol(mv1, mv2)
    T = promote_type(eltype(mv1), eltype(mv2))
    sqrt(eps(float(T)))
end

"""
    isapprox(mv1::MultiVector, mv2::MultiVector; kw...)

Check if `mv1` and `mv2` belong to the same algebra and their coefficients
are close.
"""
function Base.isapprox(mv1::MultiVector, mv2::MultiVector; 
        atol=default_atol(mv1, mv2),
        rtol=default_rtol(mv1, mv2),
)
    algebra(mv1) === algebra(mv2) || return false
    iszero(mv1) && iszero(mv2) && return true
    n1 = norm(coefficients(mv1))
    n2 = norm(coefficients(mv2))
    n12 = norm(coefficients(mv1 - mv2))
    n12 < max(atol, rtol*max(n1,n2))
end
function Base.isapprox(mv::MultiVector, x::Real; kw...)
    CA = Algebra(mv)
    isapprox(mv, MultiVector(CA,x); kw...)
end
function Base.isapprox(x::Real, mv::MultiVector; kw...)
    CA = Algebra(mv)
    isapprox(mv, MultiVector(CA,x); kw...)
end

"""
    norm_sqr(::MultiVector)

Calculates the MultiVector squared norm defined as grade(mv*reverse(mv),0)
"""
@generated function norm_sqr(mv::MultiVector{CA,T,BI,K}) where {CA,T,BI,K}
    bt = basetable(CA)
    terms = Expr[]
    for i in 1:K
        sign = baseblade_sqr(CA, bt[BI[i]])
        if !iszero(sign)
            term = :($sign*c[$i]^2)
            push!(terms, term)
        end
    end
    if isempty(terms)
        return :(zero($T))
    end
    quote
        c = coefficients(mv)
        +($(terms...))
    end
end

function baseblade_sqr(CA::Type{<:CliffordAlgebra}, x::NTuple)::Int
    prod(x, init=1) do i
        basesignature(CA,i)
    end
end

"""
    inv(::MultiVector)

Finds the inverse of the MultiVector. If no inverse exists a SingularException is thrown.
"""
@generated function inv(mv::MultiVector)
#function inv_generator(mv::Type{<:MultiVector})
    # Strategy: We can reduce the size of the problem
    # by looking only at self-adjoint multivectors.
    # Note that inv(a) == a' * inv(a' a) == a' * inv(a a').
    # Both a a' and a' a are self-adjoint with respect to reversal and
    # are therefore in a lower dimensional subspace of the algebra.
    # If b is self-adjoint, then so is inv(b) because 
    # b * inv(b) = 1 = (1)' = (inv(b) * b)' = b' * inv(b)' = b * inv(b)'
    # => inv(b) == inv(b)'
    # A multivector that is self-adjoint also has a basis expansion that
    # consists only of self-adjoint basis vectors.

    CA = algebra(mv)
    T = eltype(mv)
    BI = baseindices(mv)
    d = dimension(CA)
    basesquares = map(n -> baseproduct(CA, n, n)[2], BI)
    reversesigns = map(n -> basereverse(CA, n), BI)
    if all(iszero, basesquares)
        # singular in all cases
        return :(throw(SingularException(0)))
    end
    if length(BI) == 1
        # trivial inverse
        rc = inv(basesquares[1])
        return quote
            c = coefficients(mv)[1]
            if iszero(c)
                throw(SingularException(1))
            end
            MultiVector(Algebra(mv), $BI, ($rc * inv(c),))
        end
    end
    if all(isone, reversesigns) || all(s -> isone(-s), reversesigns)
        # We're in the (anti-)self-adjoint subspace already
        BIinv = Tuple(findall(isequal(reversesigns[1]), ntuple(i -> basereverse(CA, i), d)))
        BIone = Tuple(
            sort(
                unique(
                    skipmissing(
                        Iterators.flatten([
                            iszero(baseproduct(CA, r, c)[2]) ? missing : baseproduct(CA, r, c)[1] for r in BI,
                            c in BIinv
                        ]),
                    ),
                ),
            ),
        )
        ea = [Expr(:call, :+) for r = 1:length(BIone), c = 1:length(BIinv)]
        for (i1, bi1) in enumerate(BI), (i2, bi2) in enumerate(BIinv)
            (bi3, pc) = baseproduct(CA, bi1, bi2)
            i3 = findfirst(isequal(bi3), BIone)
            push!(ea[i3, i2].args, :($pc * c[$i1]))
        end
        if d <= 8   # this is just a guess and needs some measurements 
            ea = map(ex -> length(ex.args) == 1 ? :(zero($T)) : ex, ea)
            ma = Expr(:call, Expr(:curly, :SMatrix, length(BIone), length(BIinv)), ea[:]...)
            id = Expr(:call, :SVector, ntuple(i -> i == 1 ? one(T) : zero(T), length(BIone))...)
            return quote
                c = coefficients(mv)
                ma = $ma
                id = $id
                cinv = ma \ id
                if all(isapprox.(ma * cinv - id, zero($T); atol = 1e-8))
                    prune(MultiVector(Algebra(mv), $BIinv, Tuple(cinv)))
                else
                    throw(SingularException(2))
                end
            end
        else
            iT = typeof(inv(one(T)))
            nzea = findall(ex -> length(ex.args) > 1, ea)
            rows = map( ci -> ci[1], nzea)
            cols = map( ci -> ci[2], nzea)
            vals = map( ci -> ea[ci], nzea)
            rowsexpr = Expr(:ref, :Int, rows...)
            colsexpr = Expr(:ref, :Int, cols...)
            valsexpr = Expr(:ref, :($iT), vals... )
            rowcount = length(BIone)
            colcount = length(BIinv)
            return quote
                c = coefficients(mv)
                rows = $rowsexpr
                cols = $colsexpr
                vals = $valsexpr
                ma = sparse(rows, cols, vals, $rowcount, $colcount)
                id = sparsevec([1],[one($iT)],$rowcount)
                cinv = Array(ma) \ Array(id)
                if all(isapprox.(ma * cinv - id, zero($iT); atol = 1e-8))
                    prune(MultiVector(Algebra(mv), $BIinv, Tuple(cinv)))
                else
                    throw(SingularException(3))
                end
            end
        end
    end
    # We reduce the general case to that of self-adjoint inverses
    products = Tuple(
        (i1, i2, baseproduct(CA, b1, b2)...) for
        (i1, b1) in enumerate(BI), (i2, b2) in enumerate(BI) if
        !iszero(baseproduct(CA, b1, b2)[2]) && isone(basereverse(CA, baseproduct(CA, b1, b2)[1]))
    )
    BIsa = Tuple(sort(unique(map(p -> p[3], products))))
    coeffexpr = [Expr(:call, :+) for b in BIsa]
    for (i1, i2, bres, c) in products
        rc = basereverse(CA, BI[i2])
        push!(
            coeffexpr[findfirst(isequal(bres), BIsa)].args,
            :($(c * rc) * coefficients(mv)[$i1] * coefficients(mv)[$i2]),
        )
    end
    K = length(BIsa)
    coeffexpr = Expr(:call, Expr(:curly, :NTuple, K, T), Expr(:call, :tuple, coeffexpr...))
    quote
        @inbounds reverse(mv) * inv(MultiVector(Algebra(mv), $BIsa, $coeffexpr))
    end
end

function inv_fallback(mv::MultiVector)
    M = matrix(mv)
    e = vector(basevector(algebra(mv), :𝟏))
    c = M \ e
    bi = findall(!iszero, c)
    MultiVector(Algebra(mv), Tuple(bi), Tuple(c[bi]))
end

"""
    a / b
    (/)(a::MultiVector{CA}, b::MuliVector{CA}) where CA

Calculates the MultiVector quotient a/b by evaluating a*inv(b).
"""
(/)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = a * inv(b)
(/)(a::Real, b::MultiVector{CA}) where {CA} = a * inv(b)
(/)(a::MultiVector{CA}, b::Real) where {CA} = a * inv(b)

"""
    b \\ a
    (\\)(b::MultiVector{CA}, b::MuliVector{CA}) where CA

Calculates the MultiVector quotient a/b by evaluating inv(b)*a.
"""
(\)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = inv(a) * b
(\)(a::Real, b::MultiVector{CA}) where {CA} = inv(a) * b
(\)(a::MultiVector{CA}, b::Real) where {CA} = inv(a) * b

"""
    a ./ b

Calculates the element wise division between the MultiVectors a and b.
"""
function broadcasted(::typeof(/), a::MultiVector{CA,Ta,BI}, b::MultiVector{CA,Tb,BI})::MultiVector where {CA,Ta,Tb,BI}
    return MultiVector(CA, BI, coefficients(a) ./ coefficients(b))
end

function broadcasted(::typeof(/), a::MultiVector{CA,Ta,BIa}, b::MultiVector{CA,Tb,BIb})::MultiVector where {CA,Ta,Tb,BIa,BIb}
    @assert BIa ⊆ BIb
    v1, v2 = coefficients(a, BIb), coefficients(b, BIb)
    return MultiVector(CA, BIb, v1 ./ v2)
end

"""
    exp(::MultiVector)

Calculates the exponential function of the MultiVector defined by analytic continuation.
The generated code is automatically specialized for the sparse representation of the MultiVector.
It may take advantage of commuting base vectors and split off exponential factors. Hyperbolic, trigonometric and nilpotent solutions are recognized and handled separately.
Calling prune or grade before exp may help to find the best algorithm for the exponential evaluation.
"""

@generated function exp(mv::MultiVector)
    expr_exp(mv)
end
function expr_exp(::Type{mv}) where {mv <: MultiVector}
    # 1. Identify (optimal) subsets of base vectors so that all vectors in such a set commute with all vectors outside the set. 
    ca = algebra(mv)
    bi = baseindices(mv)
    ct = [iszero(baseproduct(ca, i, j)[2] - baseproduct(ca, j, i)[2]) for i in bi, j in bi]
    ncsets = Vector{Int}[]
    for (ki, i) in enumerate(bi)
        nc = union(collect(bi[findall(!, ct[ki, :])]), Int[i])
        si = findfirst(g -> !isempty(intersect(g, nc)), ncsets)
        if isnothing(si)
            push!(ncsets, nc)
        else
            ncsets[si] = union(nc, ncsets[si])
        end
    end
    ncsets = map(sort, ncsets)
    prodexpr = Expr(:call, :*, :(one(eltype(mv))))
    for ncset in ncsets
        if length(ncset) == 1 && ncset[1] == 1
            push!(prodexpr.args, :(exp(scalar(mv))))
        else
            baseselector = Tuple([findfirst(isequal(i), bi) for i in ncset])
            coeftuple = Expr(:tuple, map(i -> :(coefficients(mv)[$i]), baseselector)...)
            ncsetmvexpr = let
                CA = Algebra(mv) 
                BI = Tuple(ncset)
                K = length(BI)
                T = eltype(mv)
                :(MultiVector{$CA, $T, $BI, $K}($coeftuple))
            end
            basesquares = [baseproduct(ca, i, i)[2] for i in ncset]
            crosstermscancel =
                all(iszero(baseproduct(ca, i, j)[2] + baseproduct(ca, j, i)[2]) for i in ncset, j in ncset if i < j)
            if all(basesquares .== basesquares[1]) && crosstermscancel
                # We can simplify the series expansion
                if isone(basesquares[1])
                    # hyperbolic
                    push!(prodexpr.args, :(exp_hyp($ncsetmvexpr)))
                elseif isone(-basesquares[1])
                    # trigonometric
                    push!(prodexpr.args, :(exp_trig($ncsetmvexpr)))
                elseif iszero(basesquares[1])
                    # affine
                    push!(prodexpr.args, :(one(eltype(mv)) + $ncsetmvexpr))
                else
                    error("Internal error. Expected base vector to square to -1, 0 or 1")
                end
            else
                # brute force series calculation required
                push!(prodexpr.args, :(exp_fallback2($ncsetmvexpr)))
            end
        end
    end
    :(@inbounds $prodexpr)
end

function inv0(x)
    y = inv(x)
    ifelse(iszero(x), zero(y), y)
end

function exp_trig(mv::MultiVector)
    T = float(eltype(mv))
    x = sqrt(sum(coefficients(mv) .^ 2, init=zero(T)))
    s,c = sincos(x)
    return c + s * mv * inv0(x)
end

function exp_hyp(mv::MultiVector)
    T = float(eltype(mv))
    x = sqrt(sum(coefficients(mv) .^ 2, init=zero(T)))
    return cosh(x) + sinh(x) * mv * inv0(x)
end


function exp_fallback(mv::MultiVector; order::Integer = 32)
    result = one(mv)
    xpower = one(mv)
    c = inv(one(eltype(mv)))
    for n = 1:order
        c /= n
        xpower *= mv
        result += c * xpower
    end
    result
end


function exp_fallback2(mv::MultiVector; order::Integer = 16)
    s = inv(2.0^order)
    x = mv * s
    h = 1 + x * (x * (x * inv(6) + inv(2)) + 1)
    for n = 1:order
        h *= h
    end
    h
end


"""
    outermorphism(A::AbstractMatrix, mv::MultiVector)

Calculates the outermorphism f of the MultiVector defined by f(v) = Av if v is in the grade-1 subspace of the algebra.
"""
@generated function outermorphism(A::AbstractMatrix, mv::MultiVector)
    ca = Algebra(mv)
    bi = baseindices(mv)
    bt = basetable(ca)
    bv = reduce(union, bt[collect(bi)])
    if isempty(bv)
        return :(mv)
    else 
        nb = Expr(:call,:tuple)
        append!(nb.args, [:(MultiVector(ca, Tuple(A[:,$(bv[i])]))) for i in eachindex(bv)])
        s = Expr(:call,:+)
        for (k,b) in enumerate(bt[collect(bi)])
            ft = map(n -> findfirst(isequal(n), bv) ,b)
            if isempty(ft)
                push!(s.args, :(coefficients(mv)[$k]))
            elseif length(ft) == 1
                push!(s.args, :(coefficients(mv)[$k] * f[$(ft[1])]))
            else
                bp = foldl((a,b) -> Expr(:call, :∧, a, b), map(n->:(f[$n]), ft))
                push!(s.args, :(coefficients(mv)[$k] * $bp))
            end
        end
        quote
            @assert size(A,1) == size(A,2) == $(order(ca)) "The matrix A must be N×N where N is the order of the algebra."
            ca = Algebra(mv)
            f = $nb
            @inbounds $s
        end
    end
end
