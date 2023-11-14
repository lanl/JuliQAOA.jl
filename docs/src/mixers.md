# Mixers

```@index
Pages = ["mixers.md"]
```

```@docs
Mixer
MixerType
```

### Pauli X Mixers
The mixer ``H_M = \sum_i X_i`` is often referred to as the transverse field mixer, and was the original mixer introduced by [Farhi et al.](https://arxiv.org/abs/1411.4028). JuliQAOA
allows for the extension to arbitrary powers of X. 

```@docs
mixer_x(n::Int, prods::AbstractVector{Int} = [1])
```

### Grover Mixer
The Grover mixer, originally introduced [here](https://arxiv.org/abs/2006.00354), works for
both constrained and unconstrained problems and is represented by the projector onto the
feasible states: `mixer_grover` = ``|F\rangle\langle F|``.

The Grover mixer has the property that all states with the same objective value have the
same amplitude after mixing. This property can be exploited to dramatically speed up
simulating QAOA with the Grover mixer.

```@docs
mixer_grover(n::Int)
```

!!! note
    For unconstrained problems, `mixer_grover(n)` is equivalent to `mixer_x(n,0:n)/2^n`. 
    However, the underlying implementation for simulating `Grover` mixers is much faster
    than that for `X` mixers.


### XY Mixers
There are two common mixers of the form ``XX + YY``, which we refer to as the Clique and
Ring mixers. Both of these preserve Hamming weight and are therefore suited for 
constrained problems.

```@docs
mixer_clique(n::Int, k::Int; file=nothing)
mixer_ring(n::Int, k::Int; file=nothing)
```

### Novel Mixers
```@docs
mixer_general(feasible_states, m)
```