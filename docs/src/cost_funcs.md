# Cost Functions

```@index
Pages = ["cost_funcs.md"]
```

```@docs
maxcut(G::SimpleGraph, x::AbstractVector; weights=Dict())
bisection(G::SimpleGraph, x::AbstractVector; weights=Dict())
densest_ksubgraph(G::SimpleGraph, x::AbstractVector)
kvertex_cover(G::SimpleGraph, x::AbstractVector)
kSAT_instance(k::Int, n::Int, m::Int)
kSAT(instance::Array, x::Array)
sk_model(n::Int)
spin_energy(H, x)
```