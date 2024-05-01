using Base.Iterators
using SparseArrays
using LinearAlgebra
using JLD2


"""
    MixerType

Characterizes the type of mixer, used to determine the optimal statevector simulation 
algorithm.

Currently accepted types are:

- `X`: Represents any mixer composed of a symmetric sum of Pauli ``X`` operators, including 
    (but not limited to) the standard transverse field mixer.
- `Grover`: Corresponds to the Grover mixer, which projects onto the feasible subspace of 
    states.
- `General`: A placeholder for user-defined custom mixers, allowing for experimentation and 
    extension of the QAOA framework.

"""
abstract type MixerType end

struct X <: MixerType end
struct Grover <: MixerType end
struct General <: MixerType end
struct WarmStart <: MixerType end


"""
Stores metadata for the mixer

- `subtype`: Can give a string to name the mixer, e.g. 'Clique'
- `n`: the number of qubits the mixer is defined on
- `k`: `nothing` for unconstrained mixers, the Hamming weight the mixer is constrained to 
    for constrained mixers
- `N`: the total number of states the mixer operates on
- `prods`: the powers of X appearing in `X` mixers (or `nothing` for non-X mixers)
"""
struct MixerLabel{T<:MixerType}
    subtype::Union{Nothing, String}
    n::Union{Nothing, Int}
    k::Union{Nothing, Int}
    N::Union{Nothing, Int}
    prods::Union{Nothing, AbstractArray}
end

"""
    Mixer{T<:MixerType}

Stores information related to the mixing Hamiltonian ``H_M``. Fields include:

- `feasible_states`: A collection of states over which the mixer is defined and operates.
- `v`, `d`, `vinv`: Eigendecomposition of the matrix representation of ``H_M``.
- `label`: Metadata for the mixer.
- `N`: The number of feasible states, automatically determined by the length of the 
    `feasible_states` array.
- `period`: Period of the matrix (calculated from the eigenvalues `d`).
"""
struct Mixer{T<:MixerType}
    feasible_states
    N::Int
    v::Union{Nothing, AbstractMatrix, AbstractVector}
    d::Union{Nothing, AbstractVector}
    vinv::Union{Nothing, AbstractMatrix, AbstractVector}
    period::Float64
    label::MixerLabel{T}
end


"""
Convenience constructor function for `Mixer`
"""
function Mixer{T}(feasible_states, v, d::AbstractVector, vinv, label) where T
    N = length(feasible_states)
    period = get_operator_period(d)
    return Mixer{T}(feasible_states, N, v, d, vinv, period, label)
end


"""
Convenience constructor function for the `MixerLabel`.
"""
function MixerLabel{T}(;subtype=nothing,n=nothing,k=nothing,N=nothing,prods=nothing) where T
    return MixerLabel{T}(subtype, n, k, N, prods)
end


"""
Format the `MixerLabel` for output in the REPL.
"""
function print_label(l::MixerLabel{T}) where T
    ret = ""
    ret *= isnothing(l.subtype) ? "$T" : "$(l.subtype)"
    ret *= " mixer on "
    if isnothing(l.N)
        ret *= "$(l.n)-qubit states"
        ret *= isnothing(l.k) ? "" : " with Hamming weight $(l.k)"
    else
        ret *= "$(l.N) states"
    end
    if !isnothing(l.prods)
        if l.prods == [1]
            ret *= ""
        else
            ret *= " with powers $(l.prods)"
        end
    end
    return ret
end

function Base.show(io::IO, m::Mixer)
    print(io, print_label(m.label))
end


"""
X-mixers can be weighted, e.g. 0.5*mixer_x(n,[2]).
"""
function Base.:*(x::Number, mixer::Mixer{X})
    new_d = mixer.d * x
    return Mixer{X}(mixer.feasible_states, nothing, new_d, nothing, mixer.label)
end


"""
X-mixers can be divided by a number e.g. mixer_x(n,[2])/2
"""
function Base.:/(mixer::Mixer{X}, x::Number)
    new_d = mixer.d / x
    return Mixer{X}(mixer.feasible_states, nothing, new_d, nothing, mixer.label)
end


"""
X-mixers can be added, e.g. mixer_x(n,[1])+mixer_x(n,[2])
"""
function Base.:+(mixer1::Mixer{X}, mixer2::Mixer{X})
    if mixer1.N != mixer2.N
        throw(DimensionMismatch("a has length $(mixer1.N), b has length $(mixer2.N)"))
    end
    new_d = mixer1.d + mixer2.d
    new_label = combine_labels(mixer1.label, mixer2.label)
    return Mixer{X}(mixer1.feasible_states, nothing, new_d, nothing, new_label)
end


"""
Generate a new `MixerLabel` when adding X-mixers.
"""
function combine_labels(l1::MixerLabel{X}, l2::MixerLabel{X})
    if l1.n != l2.n
        throw(DimensionMismatch("l1 is for $(l1.n) qubits, l2 is for $(l2.n) qubits"))
    end
    new_prods = union(l1.prods, l2.prods)
    return MixerLabel{X}(n=l1.n,prods=new_prods)
end


"""
    mixer_x(n, prods=[1])

Create a mixer composed of a symmetric sum of products of Pauli X operators. 

`prods` determines which powers of X are included, and must be a list of ``k`` distinct 
integers ``{w_1,...,w_k}``, ``0\\le w_i \\le n``. Then `mixer(n, prods)` returns a mixer of 
the form

```math
\\sum_{j=1}^k \\sum_{i_1 < i_2 <\\ldots < i_{w_j}} X_{i_1}X_{i_2}\\cdots X_{i_{w_j}}.
```
This is a slightly unwieldy definition, and is much easier to understand by looking at 
examples

```julia-repl
julia> mixer_x(6) # = sum_i X_i
X mixer on 6-qubit states

julia> mixer_x(6, [2,3]) # = sum_{i<j} X_i X_j + sum_{i<j<k} X_i X_j X_k
X mixer on 6-qubit states with powers [2, 3]
```

Note that `X` mixer objects can be added together and given arbitrary coefficients:

```julia-repl
julia> mixer_x(6,[1]) + 0.2*mixer_x(6,[2]) + mixer_x(6, [5])/π
X mixer on 6-qubit states with powers [1, 2, 5]
```

"""
function mixer_x(n::Int, prods::AbstractVector{Int} = [1])
    d = zeros(Int, 2^n)
    for p in prods
        vals = Dict()
        for w in 0:n
            vals[w] = sum([(-1)^k*binomial(n-w,p-k)*binomial(w,k) for k in 0:w])
        end
        d += [vals[sum(s)] for s in states(n)]
    end
    label = MixerLabel{X}(n=n,prods=prods)
    return Mixer{X}(states(n), nothing, d, nothing, label)
end


"""
    mixer_grover(n)

Create the Grover mixer over all ``n``-qubit states.

    mixer_grover(n,k)

Create the Grover mixer specifically for [`dicke_states(n,k)`](@ref).
"""
function mixer_grover(n::Int)
    N = 2^n
    label = MixerLabel{Grover}(n=n)
    return Mixer{Grover}(states(n), N, nothing, nothing, nothing, 2π, label)
end

function mixer_grover(n::Int, k::Int)
    label = MixerLabel{Grover}(n=n,k=k)
    N = binomial(n,k)
    return Mixer{Grover}(dicke_states(n,k), N, nothing, nothing, nothing, 2π, label)
end

# add grover for general states
#function mixer_grover(n::Float64)
#    # to get a Grover mixer for non-standard dimension N,
#    # call mixer_grover(log2(N))
#    N = Int(round(2^n))
#    label = MixerLabel{Grover}(N=N)
#    return Mixer{Grover}(N, nothing, nothing, nothing, label)
#end

"""
    mixer_clique(n,k; file=nothing)

Create the Clique mixer ``\\frac{1}{2}\\sum_{i<j}X_i X_j + Y_i Y_j`` specifically for
[`dicke_states(n,k)`](@ref).

Works by calculating the full ``2^n \\times 2^n`` Clique mixer and then deleting rows and
columns corresponding to states with Hamming weight ``\\ne k``. The eigendecomposition of
the resulting ``\\binom{n}{k} \\times \\binom{n}{k}`` matrix is then calculated and stored.

Creating these mixers can take some time for larger `n` and `k`, in this case a 
saving+loading location `file` can be provided. If `file` already exists, `mixer_clique` 
will simply load the existing mixer in that location. If `file` does not exist, 
`mixer_clique` will calculate the mixer and store the results at `file`.

JuliQAOA uses [JLD2](https://github.com/JuliaIO/JLD2.jl) for file storage, so it is
recommended (though not necessary) that filenames passed to `file` end in `.jld2`.
"""
function mixer_clique(n::Int, k::Int; file=nothing)
    if isnothing(file) || !isfile(file)
        H = spzeros(2^n, 2^n)
        Is = [sparse(I, 2^(i-1), 2^(i-1)) for i in 1:n]
        for i in 1:n
            for j in i+1:n
                H += 0.5*kron(Is[i], _XMAT, Is[j-i], _XMAT, Is[n-j+1])
                H += 0.5*kron(Is[i], _YMAT, Is[j-i], _YMAT, Is[n-j+1])
            end
        end
        target_locs = sort([parse(Int, join(x), base=2)+1 for x in dicke_states(n,k)])
        Hreduced = collect(real.(H[target_locs, target_locs]))
        v, d, vinv = eigendecomposition(Hreduced)
        label = MixerLabel{General}(subtype="Clique",n=n,k=k)
        mixer = Mixer{General}(dicke_states(n,k), v, d, vinv, label)
        if isa(file, String)
            save_mixer(mixer, file)
        end
        return mixer
    else
        return load_mixer(file)
    end
end

"""
    mixer_ring(n,k; file=nothing)

Create the Ring mixer ``\\frac{1}{2}\\sum_i X_i X_{i+1} + Y_i Y_{i+1}`` specifically for
[`dicke_states(n,k)`](@ref).

See [`mixer_clique`](@ref) for more information and a description of the optional `file`
argument.
"""
function mixer_ring(n::Int, k::Int; file=nothing)
    if isnothing(file) || !isfile(file)
        H = spzeros(2^n, 2^n)
        Is = [sparse(I, 2^(i-1), 2^(i-1)) for i in 1:n]
        for i in 1:n-1
            H += 0.5*kron(Is[i], _XMAT, _XMAT, Is[n-i]) 
            H += 0.5*kron(Is[i], _YMAT, _YMAT, Is[n-i])
        end
        target_locs = sort([parse(Int, join(x), base=2)+1 for x in dicke_states(n,k)])
        Hreduced = collect(real.(H[target_locs, target_locs]))
        v, d, vinv = eigendecomposition(Hreduced)
        label = MixerLabel{General}(subtype="Ring",n=n,k=k)
        mixer = Mixer{General}(dicke_states(n,k), v, d, vinv, label)
        if isa(file, String)
            save_mixer(mixer, file)
        end
        return mixer
    else
        return load_mixer(file)
    end
end

"""
    mixer_warmstart(Rs)

Create a `Mixer` object from a list of single-qubit rotations `Rs`.

For an `n`-qubit QAOA, `Rs` must be of length `n`.
"""
function mixer_warmstart(Rs)
    n = length(Rs)
    Us = Vector{AbstractMatrix}()
    d = zeros(Int, 2^n)
    for i in 1:n
        v, u = eigen(Rs[i])
        push!(Us, u)
        d += kron([ones(2^(i-1)), v, ones(2^(n-i))]...)
    end
    label = MixerLabel{WarmStart}(n=n)
    return Mixer{WarmStart}(states(n), Us, d, [u' for u in Us], label)
end

"""
    mixer_general(feasible_states, m; file=nothing, name="")

Create a `Mixer` object from a list of `feasible_states` and a matrix `m`. 

Optional arguments include:
- `file`: Location for saving the mixer
- `name`: A name for the mixer
"""
function mixer_general(feasible_states, m; file=nothing, name="")
    if size(m, 1) != size(m, 2)
        throw(DimensionMismatch("Mixing matrix is not square"))
    end
    if isnothing(file) || !isfile(file)
        v, d, vinv = eigendecomposition(m)
        if isempty(name)
            label = MixerLabel{General}(N=size(m,1))
        else
            label = MixerLabel{General}(subtype=name, N=size(m,1))
        end
        mixer = Mixer{General}(feasible_states, v, d, vinv, label)
        if !isnothing(file)
            save_mixer(mixer, file)
        end
        return mixer
    else
        return load_mixer(file)
    end
end

function eigendecomposition(m::Matrix)
    d, v = eigen(m)
    if all(abs.(m - m') .< 1e-10) # check if Hermitian and ignore floating point errors
        return v, d, v'
    else
        return v, d, inv(v)
    end
end

function load_mixer(file::String)
    @load file mixer
    return mixer
end

function save_mixer(mixer::Mixer, file::String)
    @save file mixer
end