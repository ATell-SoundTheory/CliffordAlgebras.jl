# Algebra generator for CliffordAlgebras.jl

import Base.show
import Combinatorics.levicivita
import PrettyTables.pretty_table


"""
    basevectorproduct(Npos::Integer, Nneg::Integer, Nzero::Integer, base::Tuple, kl::Integer, kr::Integer)

Finds the product of the kl-th base vector with the kr-th base vector. Returns a base vector index and a scale factor.
"""
function basevectorproduct(Npos::Integer, Nneg::Integer, base::Tuple, kl::Integer, kr::Integer)
    K = length(base)
    @assert kr <= K && kl <= K
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
    (bi, coeff)
end


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
        M[kl, kr] = basevectorproduct(Npos, Nneg, base, kl, kr)
    end
    ntuple(row -> ntuple(col -> M[row, col], K), K)
end

"""
    CliffordAlgebra(Npos::Integer, Nneg::Integer, Nzero::Integer, S::NTuple(N,Symbol))

Singleton instance of the type CliffordAlgebra that describes a geometric algebra with the signature (Npos,Nneg,Nzero), base symbols S.
The base symbols are in order of the signature.
"""
struct CliffordAlgebra{Np,Nn,Nz,S,BT}
    function CliffordAlgebra(
        Npos::Integer,
        Nneg::Integer,
        Nzero::Integer,
        BaseSymbols::NTuple{N,Symbol},
    ) where {N}
        Npos = Int(Npos)
        Nneg = Int(Nneg)
        Nzero = Int(Nzero)
        @assert Npos >= 0 && Nneg >= 0 && Nzero >=0 "Algebra signature must be non-negative."
        @assert Npos + Nneg + Nzero == N "Base symbol count must match signature."
        BT = adaptbasefordual(enumeratebase(Vector, Int(N)))
        new{Npos,Nneg,Nzero,BaseSymbols,BT}()
    end
end

"""
    CliffordAlgebra(N::Integer)

Generates a geometric algebra with signature (N,0,0).
"""
function CliffordAlgebra(N::Integer)
    @assert N>=0 "Algbra signature must be non-negative."
    CliffordAlgebra(N, 0, 0, ntuple(i -> Symbol(:e, i), N))
end

"""
    CliffordAlgebra(Npos::Integer, Nneg::Integer)

Generates a geometric algebra with signature (Npos,Nneg,0).
"""
function CliffordAlgebra(Npos::Integer, Nneg::Integer)
    @assert Npos >= 0 && Nneg >= 0 "Algbra signature must be non-negative."
    CliffordAlgebra(Npos, Nneg, 0, ntuple(i -> Symbol(:e, i), Npos + Nneg))
end

"""
    CliffordAlgebra(Npos::Integer, Nneg::Integer, Nzero::Integer)

Generates a geometric algebra with signature (Npos,Nneg,Nzero).
"""
function CliffordAlgebra(Npos::Integer, Nneg::Integer, Nzero::Integer)
    @assert Npos >= 0 && Nneg >= 0 && Nzero >= 0 "Algbra signature must be non-negative."
    CliffordAlgebra(Npos, Nneg, Nzero, ntuple(i -> Symbol(:e, i), Npos + Nneg + Nzero))
end

"""
    CliffordAlgebra(a::Symbol)

Generates a predefined algebra from a identifier. Known algebras are
    - :Hyperbolic or :Hyper
    - :Complex or :ℂ
    - :Dual or :Grassmann
    - :Grassmann2D or :G2
    - :Grassmann3D or :G3
    - :Quaternions or :ℍ
    - :Cl2 and :Cl3
    - :Spacetime
    - :PGA2D or :Projective2D or :Plane2D
    - :PGA3D or :Projective3D or :Plane3D
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
    elseif a in (:Grassmann2D, :G2)
        return CliffordAlgebra(0, 0, 2, (:ε₁, :ε₂))
    elseif a in (:Grassmann3D, :G3)
        return CliffordAlgebra(0, 0, 3, (:ε₁, :ε₂, :ε₃))
    elseif a in (:Quaternions, :ℍ)
        return CliffordAlgebra(0, 2, 0, (:i, :j))
    elseif a in (:Cl2,)
        return CliffordAlgebra(2)
    elseif a in (:Cl3,)
        return CliffordAlgebra(3)
    elseif a in (:Spacetime, :STA)
        return CliffordAlgebra(1, 3, 0, (:t, :x, :y, :z))
    elseif a in (:PGA2D, :Projective2D, :Plane2D)
        return CliffordAlgebra(2, 0, 1, (:e1, :e2, :e0))
    elseif a in (:PGA3D, :Projective3D, :Plane3D)
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
        return CliffordAlgebra(4, 8, 0, (:t₁, :t₂, :e₊₁, :e₊₂, :x₁, :x₂, :y₁, :y₂, :z₁, :z₂, :e₋₁, :e₋₂))
    elseif a in (:QCGA, :QuadricConformal)
        return CliffordAlgebra(9, 6)
    else
        throw(ArgumentError("Unknown algebra."))
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
multtable(::Type{<:CliffordAlgebra{Np,Nn,Nz,S,BT}}) where {Np,Nn,Nz,S,BT} = multiplicationtable(Np, Nn, Nz, BT)
multtable(ca::CliffordAlgebra) = multtable(typeof(ca))

"""
    baseproduct(::CliffordAlgebra, nleft::Integer, nright::Integer)
    baseproduct(::Type{<:CliffordAlgebra}, nleft::Integer, nright::Integer)

Returns a tuple holding the basis index and the scalar coefficient of the geometric product of the two basis vectors with indices nleft and nright.
"""
baseproduct(::Type{CliffordAlgebra{Np,Nn,Nz,S,BT}}, nleft::Integer, nright::Integer) where {Np,Nn,Nz,S,BT} = basevectorproduct(Np, Nn, BT, nleft, nright)
baseproduct(ca::CliffordAlgebra, nleft::Integer, nright::Integer) = baseproduct(typeof(ca), nleft, nright)

"""
    order(::CliffordAlgebra)
    order(::Type{<:CliffordAlgebra})

Returns the order of the algebra. The order is the sum of the signature and the dimension of the underlying 1-vector space and the maximum grade for multivectors.
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
function basesignature(CA::Type{<:CliffordAlgebra{Np,Nn,Nz}}, n::Integer) where {Np,Nn,Nz}
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
character(CA::Type{<:CliffordAlgebra}) = baseproduct(CA, dimension(CA), dimension(CA))[2]
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


"""
    caleytable(io::IO, ca::CliffordAlgebra)
    caleytable(io::IO, CA::Type{<:CliffordAlgebra})

Generates a Cayley table view of the algebra.
"""
cayleytable(io::IO, ca::CliffordAlgebra) = cayleytable(io, typeof(ca))

function cayleytable(io::IO, CA::Type{<:CliffordAlgebra})
    dim = dimension(CA)
    bs = map(s -> s == Symbol() ? Symbol("1") : s, ntuple(k -> basesymbol(CA, k), dim))
    mt = multtable(CA)
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
    lines = vcat(0, cumsum([binomial(order(CA), i) for i = 0:order(CA)]))
    pretty_table(
        io,
        table;
        show_row_number = false,
        crop = :none,
        show_header = false,
        alignment = :c,
        vlines = lines,
        hlines = lines,
    )
end


"""
    signaturetable(io::IO, ca::CliffordAlgebra)
    signaturetable(io::IO, CA::Type{<:CliffordAlgebra})

Prints the 1-vector basis symbols and their squares.
"""
signaturetable(io::IO, ca::CliffordAlgebra) = signaturetable(io, typeof(ca))

function signaturetable(io::IO, CA::Type{<:CliffordAlgebra})
    bs = ntuple( k -> basesymbol(CA, k+1), order(CA))
    table = Matrix{String}(undef, order(CA), 2)
    for (n,s) = enumerate(bs)
        sig = basesignature(CA,n)
        table[n,1] = string(s)
        table[n,2] = sig == +1 ? "+1" : sig == -1 ? "-1" : "0"
    end
    pretty_table(
        io,
        table;
        show_row_number = false,
        crop = :none,
        show_header = false,
        alignment = :c
    )
end

function show(io::IO, ca::CliffordAlgebra)
    (Np, Nn, Nz) = signature(ca)
    print(io, "Cl(", Np, ",", Nn, ",", Nz, ")")
end
