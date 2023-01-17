# Utility functions for CliffordAlgebras.jl

"""
    nextselector(s::NTuple{N,<:Integer},M) where N

Helper function for enumerating base vectors.
"""
function nextselector(s::NTuple{N,<:Integer}, M) where {N}
    @assert issorted(s)
    @assert all(s .<= M) && all(s .>= 1)
    for n = 1:(N-1)
        if s[n] + 1 < s[n+1]
            return ntuple(i -> i < n ? i : i == n ? s[i] + 1 : s[i], N)
        end
    end
    if s[N] < M
        return ntuple(i -> i < N ? i : i == N ? s[i] + 1 : s[i], N)
    else
        return nothing
    end
end

"""
    enumeratebase(N::Integer)

Generates the basis vector product combinations for the Clifford algebra of an N dimensional vector space.
"""
function enumeratebase(N::Integer)::Tuple
    Tuple(enumeratebase(Vector, N))
end
function enumeratebase(::Type{Vector}, N::Integer)::Vector
    result = Tuple[]
    push!(result, ())
    for n = 1:N
        s = ntuple(identity, n)
        while true
            push!(result, s)
            s = nextselector(s, N)
            !isnothing(s) || break
        end
    end
    @assert length(result) == 2^N
    return result
end


"""
    adaptbasefordual(base)

Shuffles the base vector products so that the Poincaré dual of every basis vector is another basis vector.
"""
function adaptbasefordual(base)::Tuple
    N = length(base)
    Nh = N ÷ 2
    Tuple(
        n > Nh && length(base[n]) >= 2 && levicivita([base[N-n+1]..., base[n]...]) < 0 ?
        Tuple(
            k == length(base[n]) ? base[n][k-1] :
            k + 1 == length(base[n]) ? base[n][k+1] : base[n][k] for k = 1:length(base[n])
        ) : base[n] for n = 1:N
    )
end

"""
    findpermutation(a,b)

Find the permutation that takes the a to b.
"""
function findpermutation(a, b)
    @assert allunique(a) "Permutations can only be found in sequences of unique elements."
    @assert length(a) == length(b) "Permutations can only be found between sequences of identical length."
    @assert length(intersect(a, b)) == length(a) "Permutations can only be found between sequences with the same elements."
    [findfirst(isequal(a[n]), b) for n in eachindex(a)]
end

