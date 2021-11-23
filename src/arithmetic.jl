# MultiVector arithmetic for CliffordAlgebras.jl

import Base.+, Base.-, Base.*, Base./, Base.\
import Base.inv, Base.adjoint, Base.exp
import LinearAlgebra.norm, LinearAlgebra.norm_sqr
import StaticArrays.SVector, StaticArrays.SMatrix


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
        return :(MultiVector(algebra(a), zero($T)))
    end
    coeffs = (
        iszero(acc_a[n]) ? :(b.c[$(acc_b[n])]) :
        iszero(acc_b[n]) ? :(a.c[$(acc_a[n])]) : :(a.c[$(acc_a[n])] + b.c[$(acc_b[n])]) for
        n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(MultiVector(algebra(a), $BI, $coeffsexpr))
end

(+)(a::Real, b::MultiVector{CA}) where {CA} = MultiVector(CA, a) + b
(+)(a::MultiVector{CA}, b::Real) where {CA} = a + MultiVector(CA, b)

(-)(a::MultiVector) = -1 * a
(-)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = a + (-b)
(-)(a::Real, b::MultiVector) = a + (-b)
(-)(a::MultiVector, b::Real) = a + (-b)


geometricfilter(leftgrade, rightgrade, productgrade) = true
exteriorfilter(leftgrade, rightgrade, productgrade) = productgrade == leftgrade + rightgrade
leftcontractionfilter(leftgrade, rightgrade, productgrade) =
    productgrade == rightgrade - leftgrade
rightcontractionfilter(leftgrade, rightgrade, productgrade) =
    productgrade == leftgrade - rightgrade
fatdotfilter(leftgrade, rightgrade, productgrade) =
    productgrade == abs(leftgrade - rightgrade)
scalarfilter(leftgrade, rightgrade, productgrade) = iszero(productgrade)

@generated function filteredprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
    ::Val{Filter},
) where {CA,Ta,Tb,Filter}
    d = dimension(CA)
    mt = multtable(CA)
    acc = [NTuple{3,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)), (kb, base_b) in enumerate(baseindices(b))
        baseidx, coeff = mt[base_a][base_b]
        if !iszero(coeff) &&
           Filter(basegrade(CA, base_a), basegrade(CA, base_b), basegrade(CA, baseidx))
            push!(acc[baseidx], (ka, kb, coeff))
        end
    end
    BI = Tuple(findall(!isempty, acc))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (:($(tup[3]) * a.c[$(tup[1])] * b.c[$(tup[2])]) for tup in acc[n])...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(MultiVector(algebra(a), $BI, $coeffsexpr))
end

"""
    a * b

Calculates the geometric product of two MultiVectors a and b.
"""
(*)(a::MultiVector, b::MultiVector) where {CA} = filteredprod(a, b, Val(geometricfilter))
(*)(a::Real, b::MultiVector) = MultiVector(algebra(b), baseindices(b), a .* coefficients(b))
(*)(a::MultiVector, b::Real) = b * a

"""
    a ∧ b

Calculates the wedge product between two MultiVectors a and b.
"""
(∧)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} =
    filteredprod(a, b, Val(exteriorfilter))

"""
    a ⋅ b

Calculates the "fat dot" product between the MultiVectors a and b.
"""
(⋅)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} =
    filteredprod(a, b, Val(fatdotfilter))

"""
    a ⨼ b

Calculates the left contraction of the MultiVectors a and b.
"""
(⨼)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} =
    filteredprod(a, b, Val(leftcontractionfilter))

"""
    a ⨽ b

Calculates the right contraction of the MultiVectors a and b.
"""
(⨽)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} =
    filteredprod(a, b, Val(rightcontractionfilter))

"""
    a ⋆ b

Calculates the scalar product of the MultiVectors a and b.
"""
(⋆)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} =
    filteredprod(a, b, Val(scalarfilter))

"""
    a ∨ b

Calculates the vee product of the MultiVectors a and b.
"""
(∨)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = dual(dual(a) ∧ dual(b))

@generated function symprod(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
    ::Val{P},
) where {CA,Ta,Tb,P}
    d = dimension(CA)
    mt = multtable(CA)
    acc = [NTuple{3,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)), (kb, base_b) in enumerate(baseindices(b))
        baseidx1, coeff1 = mt[base_a][base_b]
        baseidx2, coeff2 = mt[base_b][base_a]
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
        return :(MultiVector(algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (:($(tup[3]) * a.c[$(tup[1])] * b.c[$(tup[2])]) for tup in acc[n])...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(MultiVector(algebra(a), $BI, $coeffsexpr))
end

"""
    a ×₋ b

Calculates the commutator ab-ba of two MultiVectors a and b.
"""
(×₋)(a::MultiVector{CA,Ta}, b::MultiVector{CA,Tb}) where {CA,Ta,Tb} =
    inv(promote_type(Ta, Tb)(2)) * symprod(a, b, Val(-1))

(×₋)(a::Ta, b::MultiVector{CA,Tb}) where {CA,Ta<:Real,Tb} = zero(promote_type(Ta, Tb))
(×₋)(a::MultiVector, b::Real) = b ×₋ a
(×₋)(a::Ta, b::Tb) where {Ta<:Real,Tb<:Real} = zero(promote_type(Ta, Tb))

"""
    a ×₊ b

Calculates the anti-commutator ab+ba of two MultiVectors a and b.
"""
(×₊)(a::MultiVector{CA,Ta}, b::MultiVector{CA,Tb}) where {CA,Ta,Tb} =
    inv(promote_type(Ta, Tb)(2)) * symprod(a, b, Val(+1))
(×₊)(a::Real, b::MultiVector) = a * b
(×₊)(a::MultiVector, b::Real) = a * b
(×₊)(a::Real, b::Real) = a * b


@generated function sandwichproduct(
    a::MultiVector{CA,Ta},
    b::MultiVector{CA,Tb},
) where {CA,Ta,Tb}
    # a * b * reverse(a)
    d = dimension(CA)
    mt = multtable(CA)
    acc = [NTuple{4,Int}[] for n = 1:d]
    for (ka, base_a) in enumerate(baseindices(a)),
        (kb, base_b) in enumerate(baseindices(b)),
        (ka_rev, base_a_rev) in enumerate(baseindices(a))

        leftbaseidx, leftcoeff = mt[base_a][base_b]
        baseidx, rightcoeff = mt[leftbaseidx][base_a_rev]
        coeff = leftcoeff * rightcoeff * basereverse(CA, base_a_rev)
        if !iszero(coeff)
            push!(acc[baseidx], (ka, kb, ka_rev, coeff))
        end
    end
    BI = Tuple(findall(!isempty, acc))
    K = length(BI)
    T = promote_type(Ta, Tb)
    if iszero(K)
        return :(MultiVector(algebra(a), zero($T)))
    end
    coeffs = (
        Expr(
            :call,
            :+,
            (
                :($(tup[4]) * a.c[$(tup[1])] * b.c[$(tup[2])] * a.c[$(tup[3])]) for
                tup in acc[n]
            )...,
        ) for n in BI
    )
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(MultiVector(algebra(a), $BI, $coeffsexpr))
end

"""
    a ≀ b

Calculates the sandwich product a*b*reverse(a) for two MultiVectors a and b.
"""
(≀)(a::MultiVector{CA}, b::MultiVector{CA}) where {CA} = sandwichproduct(a, b)

"""
    polarize(mv::MultiVector)
    mv'

Calculates the polarization of the MultiVector, i.e. mv * pseudoscalar.
"""
polarize(mv::MultiVector{CA}) where {CA} = mv * pseudoscalar(CA)

adjoint(mv::MultiVector) = polarize(mv)


"""
    Λᵏ(::MultiVector, ::Val{k}) where k
    Λᵏ(::MultiVector, k::Integer) 

Projects the MultiVector onto k-vectors. Similar to grade(mv,k), but uses
@generated code and compile time optimizations.
"""
@generated function (Λᵏ)(mv::MultiVector, ::Val{k}) where {k}
    CA = algebra(mv)
    indexbounds =
        1 .+ (0, cumsum(ntuple(i -> binomial(order(CA), i - 1), 1 + order(CA)))...)
    @assert 0 < k + 1 < length(indexbounds) "k out of bounds."
    bi_low = indexbounds[k+1]
    bi_high = indexbounds[k+2] - 1
    s = findall(i -> i >= bi_low && i <= bi_high, baseindices(mv))
    BI = baseindices(mv)[s]
    K = length(BI)
    T = eltype(mv)
    coeffs = (:(mv.c[$i]) for i in s)
    if length(coeffs) == 0
        return :(MultiVector(algebra(a), zero($T)))
    end
    coeffsexpr = Expr(:call, :(NTuple{$K,$T}), Expr(:call, :tuple, coeffs...))
    :(MultiVector(algebra(mv), $BI, $coeffsexpr))
end

(Λᵏ)(mv::MultiVector, k::Integer) = Λᵏ(mv, Val(k))

"""
    norm(::MultiVector)

Calculates the MultiVector norm defined as sqrt(grade(mv*reverse(mv),0))
"""
norm(mv::MultiVector) = sqrt(scalar(Λᵏ(mv * reverse(mv), 0)))

"""
    norm_sqr(::MultiVector)

Calculates the MultiVector squared norm defined as grade(mv*reverse(mv),0)
"""
norm_sqr(mv::MultiVector) = scalar(Λᵏ(mv * reverse(mv), 0))


"""
    inv(::MultiVector)

Finds the inverse of the MultiVector. If no inverse exists a SingularException is thrown.
"""
@generated function inv(mv::MultiVector)
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
    mt = multtable(CA)
    basesquares = map(n -> mt[n][n][2], BI)
    reversesigns = map(n -> basereverse(CA, n), BI)
    if all(iszero, basesquares)
        # singular in all cases
        return :(throw(SingularException))
    end
    if length(BI) == 1
        # trivial inverse
        rc = inv(basesquares[1])
        return quote
            c = coefficients(mv)[1]
            if iszero(c)
                throw(SingularException)
            end
            MultiVector(algebra(mv), $BI, ($rc * inv(c),))
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
                            iszero(mt[r][c][2]) ? missing : mt[r][c][1] for r in BI,
                            c in BIinv
                        ]),
                    ),
                ),
            ),
        )
        ea = [Expr(:call, :+) for r = 1:length(BIone), c = 1:length(BIinv)]
        for (i1, bi1) in enumerate(BI), (i2, bi2) in enumerate(BIinv)
            (bi3, pc) = mt[bi1][bi2]
            i3 = findfirst(isequal(bi3), BIone)
            push!(ea[i3, i2].args, :($pc * c[$i1]))
        end
        ea = map(ex -> length(ex.args) == 1 ? :(zero($T)) : ex, ea)
        ma = Expr(:call, Expr(:curly, :SMatrix, length(BIone), length(BIinv)), ea[:]...)
        id = Expr(:call, :SVector, ntuple(i -> i == 1 ? one(T) : zero(T), length(BIone))...)
        return quote
            c = coefficients(mv)
            cinv = $ma \ $id
            if all(isapprox.($ma * cinv - $id, zero($T); atol = 1e-8))
                MultiVector(algebra(mv), $BIinv, Tuple(cinv))
            else
                throw(SingularException)
            end
        end
    end
    # We reduce the general case to that of self-adjoint inverses
    products = Tuple(
        (i1, i2, mt[b1][b2][1], mt[b1][b2][2]) for
        (i1, b1) in enumerate(BI), (i2, b2) in enumerate(BI) if
        !iszero(mt[b1][b2][2]) && isone(basereverse(CA, mt[b1][b2][1]))
    )
    BIsa = Tuple(sort(unique(map(p -> p[3], products))))
    coeffexpr = [Expr(:call, :+) for b in BIsa]
    for (i1, i2, bres, c) in products
        rc = basereverse(CA, BI[i2])
        push!(
            coeffexpr[findfirst(isequal(bres), BIsa)].args,
            :($(c * rc) * mv.c[$i1] * mv.c[$i2]),
        )
    end
    K = length(BIsa)
    coeffexpr = Expr(:call, Expr(:curly, :NTuple, K, T), Expr(:call, :tuple, coeffexpr...))
    quote
        reverse(mv) * inv(MultiVector(algebra(mv), $BIsa, $coeffexpr))
    end
end

function inv_fallback(mv::MultiVector)
    M = matrix(mv)
    e = vector(basevector(algebra(mv), :𝟏))
    c = M \ e
    bi = findall(!iszero, c)
    MultiVector(algebra(mv), Tuple(bi), Tuple(c[bi]))
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
    exp(::MultiVector)

Calculates the exponential function of the MultiVector defined by analytic continuation.
The generated code is automaticall specialized for the sparse representation of the MultiVector.
It may take advantage of commuting base vectors and split off exponential factors. Hyperbolic, trigonometric and nilpotent solutions are recognized and handled separately.
Calling prune or grade before exp may help to find the best algorithm for the exponential evaluation.
"""
@generated function exp(mv::MultiVector)
    # 1. Identify (optimal) subsets of base vectors so that all vectors in such a set commute with all vectors outside the set. 
    ca = algebra(mv)
    bi = baseindices(mv)
    mt = multtable(ca)
    ct = [iszero(mt[i][j][2] - mt[j][i][2]) for i in bi, j in bi]
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
            coeftuple = Expr(:call, :tuple, map(i -> :(mv.c[$i]), baseselector)...)
            ncsetmvexpr = :(MultiVector(algebra(mv), $(Tuple(ncset)), $coeftuple))
            basesquares = [mt[i][i][2] for i in ncset]
            crosstermscancel =
                all(iszero(mt[i][j][2] + mt[j][i][2]) for i in ncset, j in ncset if i < j)
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
    prodexpr
end


function exp_trig(mv::MultiVector)
    if iszero(mv)
        return one(mv)
    else
        x = sqrt(sum(coefficients(mv) .^ 2))
        return cos(x) + sin(x) * mv * inv(x)
    end
end


function exp_hyp(mv::MultiVector)
    if iszero(mv)
        return one(mv)
    else
        x = sqrt(sum(coefficients(mv) .^ 2))
        return cosh(x) + sinh(x) * mv * inv(x)
    end
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
    ca = algebra(mv)
    bi = baseindices(mv)
    bt = basetable(ca)
    bv = reduce(union, bt[collect(bi)])
    if isempty(bv)
        return :(mv)
    else 
        nb = Expr(:call,:tuple)
        append!(nb.args, [:(MultiVector(ca, Tuple(A[:,$(bv[i])]))) for i = 1:length(bv)])
        s = Expr(:call,:+)
        for (k,b) in enumerate(bt[collect(bi)])
            ft = map(n -> findfirst(isequal(n), bv) ,b)
            if isempty(ft)
                push!(s.args, :(mv.c[$k]))
            elseif length(ft) == 1
                push!(s.args, :(mv.c[$k] * f[$(ft[1])]))
            else
                bp = foldl((a,b) -> Expr(:call, :∧, a, b), map(n->:(f[$n]), ft))
                push!(s.args, :(mv.c[$k] * $bp))
            end
        end
        quote
            ca = algebra(mv)
            f = $nb
            $s
        end
    end
end