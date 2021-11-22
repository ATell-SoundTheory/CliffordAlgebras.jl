# Algebra generator for CliffordAlgebras.jl

import Base.show
import Combinatorics.levicivita
import PrettyTables.pretty_table

"""
    multiplicationstable(Npos::Integer, Nneg::Integer, Nzero::Integer, base::Tuple)

Generates a multiplication table for the Clifford algebra with the specified signature and base vectors.
"""
function multiplicationtable(Npos::Integer, Nneg::Integer, Nzero::Integer, base::Tuple)
    N = Npos + Nneg + Nzero
    K = 2^N
    @assert length(base) == K
    M = Matrix{Tuple}(undef, K, K)
    for kl = 1:K, kr = 1:K
        selector_left = base[kl]
        selector_right = base[kr]
        selector_prod = collect((selector_left..., selector_right...))
        permutation = sortperm(selector_prod)
        coeff = levicivita(permutation)
        v = []
        for b in unique(selector_prod[permutation])
            c = count(isequal(b), selector_prod)
            @assert c in (1, 2)
            if c == 1
                push!(v, b)
            else
                coeff *= (b <= Npos) ? 1 : (b <= Npos + Nneg) ? -1 : 0
            end
        end
        v = Tuple(v)
        vl = length(v)
        bi = findfirst(b -> length(b) == vl && length(intersect(b, v)) == vl, base)
        coeff *= levicivita(findpermutation(v, base[bi]))
        M[kl, kr] = (bi, coeff)
    end
    ntuple(row -> ntuple(col -> M[row, col], K), K)
end

"""
    CliffordAlgebra

CliffordAlgebra{Np,Nn,Nz,S} is a type that describes a geometric algebra with the signature (Np,Nn,Nz), base symbols S.
"""
struct CliffordAlgebra{Np,Nn,Nz,S,BT,MT}
    function CliffordAlgebra(
        Npos::Integer,
        Nneg::Integer,
        Nzero::Integer,
        BaseSymbols::NTuple{N,Symbol},
    ) where {N}
        Npos = Int(Npos)
        Nneg = Int(Nneg)
        Nzero = Int(Nzero)
        @assert Npos + Nneg + Nzero == N
        BT = adaptbasefordual(enumeratebase(Int(N)))
        MT = multiplicationtable(Npos, Nneg, Nzero, BT)
        new{Npos,Nneg,Nzero,BaseSymbols,BT,MT}()
    end
end

"""
    CliffordAlgebra(N::Integer)

Generates a geometric algebra with signature (N,0,0).
"""
function CliffordAlgebra(N::Integer)
    CliffordAlgebra(N, 0, 0, ntuple(i -> Symbol(:e, i), N))
end

"""
    CliffordAlgebra(Npos::Integer, Nneg::Integer)

Generates a geometric algebra with signature (Npos,Nneg,0).
"""
function CliffordAlgebra(Npos::Integer, Nneg::Integer)
    CliffordAlgebra(Npos, Nneg, 0, ntuple(i -> Symbol(:e, i), Npos + Nneg))
end

"""
    CliffordAlgebra(Npos::Integer, Nneg::Integer, Nzero::Integer)

Generates a geometric algebra with signature (Npos,Nneg,Nzero).
"""
function CliffordAlgebra(Npos::Integer, Nneg::Integer, Nzero::Integer)
    CliffordAlgebra(Npos, Nneg, Nzero, ntuple(i -> Symbol(:e, i), Npos + Nneg + Nzero))
end

"""
    CliffordAlgebra(a::Symbol)

Generates a predefined algebra from a identifier. Known algebras are
    - :Hyperbolic or :Hyper
    - :Complex or :ℂ
    - :Dual or :Grassmann
    - :Quaternions or :ℍ
    - :Cl2 and :Cl3
    - :Spacetime
    - :PGA2D or :Projective2D
    - :PGA3D or :Projective3D
    - :CGA2D or :Conformal2D
    - :CGA3D or :Conformal3D
    - :DCGA3D or :DoubleConformal3D
    - :TCGA3D or :TripleConformal3D
    - :DCGSTA or :DoubleConformalSpacetime
    - :QCGA or :QuadricConformal
"""
function CliffordAlgebra(a::Symbol)
    if a in (:Hyperbolic, :Hyper)
        return CliffordAlgebra(1, 0, 0, (:j,))
    elseif a in (:Complex, :ℂ)
        return CliffordAlgebra(0, 1, 0, (:i,))
    elseif a in (:Dual, :Grassmann)
        return CliffordAlgebra(0, 0, 1, (:ε,))
    elseif a in (:Quaternions, :ℍ)
        return CliffordAlgebra(0, 2, 0, (:i, :j))
    elseif a in (:Cl2,)
        return CliffordAlgebra(2)
    elseif a in (:Cl3,)
        return CliffordAlgebra(3)
    elseif a in (:Spacetime, :STA)
        return CliffordAlgebra(1, 3, 0, (:t, :x, :y, :z))
    elseif a in (:PGA2D, :Projective2D)
        return CliffordAlgebra(2, 0, 1, (:e1, :e2, :e0))
    elseif a in (:PGA3D, :Projective3D)
        return CliffordAlgebra(3, 0, 1, (:e1, :e2, :e3, :e0))
    elseif a in (:CGA2D, :Conformal2D)
        return CliffordAlgebra(3, 1, 0, (:e1, :e2, :e₊, :e₋))
    elseif a in (:CGA3D, :Conformal3D)
        return CliffordAlgebra(4, 1, 0, (:e1, :e2, :e3, :e₊, :e₋))
    elseif a in (:DCGA3D, :DoubleConformal3D)
        return CliffordAlgebra(6, 2)
    elseif a in (:TCGA3D, :TripleConformal3D)
        return CliffordAlgebra(9, 3)
    elseif a in (:DCGSTA, :DoubleConformalSpacetime)
        return CliffordAlgebra(4, 8)
    elseif a in (:QCGA, :QuadricConformal)
        return CliffordAlgebra(9, 6)
    else
        error("Unknown algebra.")
    end
end


"""
    basesymbols(::CliffordAlgebra)
    basesymbols(::Type{<:CliffordAlgebra})

Returns the 1-vector space basis symbols of the algebra.
"""
basesymbols(::Type{<:CliffordAlgebra{Np,Nn,Nz,S}}) where {Np,Nn,Nz,S} = S
basesymbols(ca::CliffordAlgebra) = basesymbols(typeof(ca))

"""
    basetable(::CliffordAlgebra)
    basetable(::Type{<:CliffordAlgebra})

Returns the internal basis table of the algebra.
"""
basetable(::Type{<:CliffordAlgebra{Np,Nn,Nz,S,BT}}) where {Np,Nn,Nz,S,BT} = BT
basetable(ca::CliffordAlgebra) = basetable(typeof(ca))

"""
    multtable(::CliffordAlgebra)
    multtable(::Type{<:CliffordAlgebra})

Returns the internal multuplication table of the geometric product of the algebra.
"""
multtable(::Type{<:CliffordAlgebra{Np,Nn,Nz,S,BT,MT}}) where {Np,Nn,Nz,S,BT,MT} = MT
multtable(ca::CliffordAlgebra) = multtable(typeof(ca))

"""
    order(::CliffordAlgebra)
    order(::Type{<:CliffordAlgebra})

Returns the order of the algebra. The order is the sum of the signature.
"""
order(::Type{<:CliffordAlgebra{Np,Nn,Nz}}) where {Np,Nn,Nz} = Np + Nn + Nz
order(ca::CliffordAlgebra) = order(typeof(ca))

"""
    signature(::CliffordAlgebra)
    signature(::Type{<:CliffordAlgebra})

Returns the signature of the algebra.
"""
signature(::Type{<:CliffordAlgebra{Np,Nn,Nz}}) where {Np,Nn,Nz} = (Np, Nn, Nz)
signature(ca::CliffordAlgebra) = signature(typeof(ca))

"""
    basesignature(::CliffordAlgebra, n::Integer)
    basesignature(::Type{<:CliffordAlgebra}, n::Integer)

Returns the signature value of the n-th basis 1-vector. The return value is +1, -1 or 0.
"""
function basesignature(::Type{<:CliffordAlgebra{Np,Nn,Nz}}, n::Integer) where {Np,Nn,Nz}
    @assert 1 <= n <= order(CA)
    if n <= Np
        return +1
    elseif n <= Np + Nn
        return -1
    else
        return 0
    end
end

basesignature(ca::CliffordAlgebra, n::Integer) = basesignature(typeof(ca), n)

"""
    dimension(::CliffordAlgebra)
    dimension(::Type{<:CliffordAlgebra})

Returns the dimension of the algebra, i.e. the number of coefficients in a general multivector.
"""
dimension(CA::Type{<:CliffordAlgebra}) = 2^order(CA)
dimension(ca::CliffordAlgebra) = dimension(typeof(ca))

"""
    character(::CliffordAlgebra)
    character(::Type{<:CliffordAlgebra})

Returns the square of the pseudoscalar of the algebra.
"""
character(CA::Type{<:CliffordAlgebra}) = multtable(CA)[dimension(CA)][dimension(CA)][2]
character(ca::CliffordAlgebra) = character(typeof(ca))

"""
    basegrade(::CliffordAlgebra, n::Integer)
    basegrade(::Type{<:CliffordAlgebra}, n::Integer)

Returns the grade of the n-th basis multivector of the algebra.
"""
basegrade(CA::Type{<:CliffordAlgebra}, n::Integer) = length(basetable(CA)[n])
basegrade(ca::CliffordAlgebra, n::Integer) = basegrade(typeof(ca), n)

"""
    basereverse(::CliffordAlgebra, n::Integer)
    basereverse(::Type{<:CliffordAlgebra}, n::Integer)

Returns the sign change of the n-th basis multivector under reversal.
"""
basereverse(CA::Type{<:CliffordAlgebra}, n::Integer) = (-1)^(basegrade(CA, n) ÷ 2)
basereverse(ca::CliffordAlgebra, n::Integer) = basereverse(typeof(ca), n)

"""
    basesymbol(::CliffordAlgebra, n::Integer)
    basesymbol(::Type{<:CliffordAlgebra}, n::Integer)

Returns the symbol used for the n-th basis multivector of the algebra.
"""
function basesymbol(CA::Type{<:CliffordAlgebra}, n::Integer)
    @assert 1 <= n <= dimension(CA)
    s = basetable(CA)[n]
    if isempty(s)
        Symbol(:𝟏)
    else
        Symbol(map(i -> basesymbols(CA)[i], s)...)
    end
end

basesymbol(ca::CliffordAlgebra, n::Integer) = basesymbol(typeof(ca), n)


function show(io::IO, ca::CliffordAlgebra)
    (Np, Nn, Nz) = signature(ca)
    println(io, "Cl(", Np, ",", Nn, ",", Nz, ")")
    dim = dimension(ca)
    bs = map(s -> s == Symbol() ? Symbol("1") : s, ntuple(k -> basesymbol(ca, k), dim))
    mt = multtable(ca)
    table = Matrix{String}(undef, dim, dim)
    for col = 1:dim
        for row = 1:dim
            baseidx, coeff = mt[row][col]
            @assert coeff in (-1, 0, 1)
            if baseidx == 1
                table[row, col] = coeff == -1 ? "-1" : coeff == 1 ? "+1" : "0"
            else
                if iszero(coeff)
                    table[row, col] = "0"
                else
                    table[row, col] = string(coeff > 0 ? :+ : :-, bs[baseidx])
                end
            end
        end
    end
    lines = vcat(0, cumsum([binomial(order(ca), i) for i = 0:order(ca)]))
    pretty_table(
        table;
        show_row_number = false,
        crop = :none,
        noheader = true,
        alignment = :c,
        vlines = lines,
        hlines = lines,
    )
end
