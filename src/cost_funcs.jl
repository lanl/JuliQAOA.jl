using Graphs
using StatsBase

"""
    maxcut(G::SimpleGraph, x; weights=Dict())

Calculate the number of edges cut by a partition `x` on graph `G`. 

`weights` can be provided as dictionary `(i,j)=>w` where `i < j`, specifying the weights for
each edge `(i,j)`.

# Example
```julia-repl
julia> using JuliQAOA, Graphs

julia> g = SimpleGraph(5)
{5, 0} undirected simple Int64 graph

julia> edges = [(1,2), (1,4), (2,3), (2,5), (3,5), (4,5)];

julia> for edge in edges
           add_edge!(g, edge[1], edge[2])
       end

julia> maxcut(g, [0, 0, 1, 1, 1])
3

julia> weights = Dict((1,2)=>4,(1,4)=>2,(2,3)=>1,(2,5)=>3,(3,5)=>3,(4,5)=>1);

julia> maxcut(g, [0, 1, 0, 0, 1]; weights=weights)
9
```
"""
function maxcut(G::SimpleGraph, x::AbstractVector; weights=Dict())
    if !all(i < j for (i, j) in keys(weights))
        throw(ArgumentError("weights must be given in terms of tuples (i,j) with i<j"))
    end
    C = 0
    for edge in Graphs.edges(G)
        u,v = src(edge), dst(edge)  
        if x[u] != x[v]
            C += get(weights, (u,v), 1)
        end
    end
    return C
end


"""
    bisection(G::SimpleGraph, x; weights=Dict())

Calculate the number of edges cut by a bisection `x` on graph `G`, and throws an error if 
`x` is not a bisection.
"""
function bisection(G::SimpleGraph, x::AbstractVector; weights=Dict())
    n = nv(G)
    if sum(x) > n/2 + 1 || sum(x) < n/2 - 1
        throw(ArgumentError("Invalid bisection size"))
    end
    return maxcut(G, x; weights=weights)
end

"""
    densest_ksubgraph(G::SimpleGraph, x)

Calculate the number of edges contained in the subgraph of `G` denoted by the vertices
marked in `x`.

# Example
```julia-repl
julia> using JuliQAOA, Graphs

julia> g = SimpleGraph(5)
{5, 0} undirected simple Int64 graph

julia> edges = [(1,2), (1,4), (2,3), (2,5), (3,5), (4,5)];

julia> for edge in edges
           add_edge!(g, edge[1], edge[2])
       end

julia> densest_ksubgraph(g, [0, 1, 1, 0, 1])
3
```
"""
function densest_ksubgraph(G::SimpleGraph, x::AbstractVector)
    C = 0
    for edge in Graphs.edges(G)
        if x[src(edge)]==x[dst(edge)]==1
            C += 1
        end
    end
    return C
end

"""
    kvertex_cover(G::SimpleGraph, x)

Calculate the number of edges touching at least one node of the subgraph of `G` denoted by 
the vertices marked in `x`.

# Example
```julia-repl
julia> using JuliQAOA, Graphs

julia> g = SimpleGraph(5)
{5, 0} undirected simple Int64 graph

julia> edges = [(1,2), (1,4), (2,3), (2,5), (3,5), (4,5)];

julia> for edge in edges
           add_edge!(g, edge[1], edge[2])
       end

julia> kvertex_cover(g, [1, 0, 1, 0, 1])
6
```
"""
function kvertex_cover(G::SimpleGraph, x::AbstractVector)
    C = 0
    for edge in Graphs.edges(G)
        if x[src(edge)]== 1 || x[dst(edge)]==1
            C += 1
        end
    end
    return C
end


"""
    kSAT_instance(k,n,m)

Generate a random `k`-SAT instance with `n` variables and `m` clauses.

Returns a list of clauses of the form e.g. `[1,-2,3]`, where the minus sign indicates
negation. To be explicit, this clause translates as "`x[1]==1` OR `x[2]==0` OR `x[3]==1`
for a variable assignment `x`. 
"""
function kSAT_instance(k::Int, n::Int, m::Int)
    instance = Vector{Vector{Int64}}()
    for _ in 1:m
        vars = sample(1:n, k, replace=false)
        clause = [(rand(Bool) ? -var : var) for var in vars] 
        push!(instance, clause)
    end
    return instance
end

"""
    kSAT(instance::Array, x)

Evaluate the number of satisfied clauses in a given kSAT `instance` for a the variable
assignment `x`. See [`kSAT_instance`](@ref) for a description of the correct format for
`instance`s.
"""
function kSAT(instance::Array, x::Array)
    C = 0
    for clause in instance
        if any(lit -> (lit > 0 && x[lit] == 1) || (lit < 0 && x[-lit] == 0), clause)
            C += 1
        end
    end
    return C
end

"""
    SK_model(n)

Generate a random Sherrington-Kirkpatrick model on `n` qubits. 

The Sherrington-Kirkpatrick model is defined by a Hamiltonian of the form

```math
H = \\frac{1}{N} \\sum_{i<j} J_{ij} Z_i Z_j
```

where the ``J_{ij}`` are sampled from a normal distribution with mean 0 and standard 
deviation 1. 

Returns a `Dict` of the form `(i,j)=>Jij`.
"""
function sk_model(n::Int)
    H = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:n
        for j in (i+1):n
            H[(i, j)] = randn()
        end
    end
    return H
end

"""
    spin_energy(H,x)

Calculate the energy of a system with Hamiltonian `H` and binary assignment `x`. 

The Hamiltonian should be provided as a `Dict` of the coupling terms. That is, the 
Hamiltonian 

```math
H = c_1 Z_1 + c_{2,3} Z_2 Z_3 + c_{4,5,6} Z_4 Z_5 Z_6 + \\ldots
```

would be input as `Dict((1,)=>c1, (2,3)=>c23, (4,5,6)=>c456, ...)`.

The binary assignment `x` is converted to a spin assignment via the standard `0 = ↑ = 1`,
`1 = ↓ = -1`.

# Example
```julia-repl
julia> using JuliQAOA

julia> H = Dict((1,)=>-1.2, (2,)=>2.7, (3,)=>1.4, (1,2)=>2.2, (1,3)=>-5.9, (1,2,3)=>-4.3);

julia> spin_energy(H, [1,0,1])
-9.9
```
"""
function spin_energy(H, x)
    σ = (-1) .^ x 
    energy = 0.0
    for (key, value) in H
        energy += prod(σ[collect(key)])*value
    end
    return energy
end