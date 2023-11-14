# Angle Finding

JuliQAOA uses [Enzyme.jl](https://enzyme.mit.edu/julia/stable/) to enable 
[automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) (also 
referred to as "AD" or "autodiff"). In short, AD allows for calculating the gradient of a
function `f(x1, x2,...)` in the same time as a single evaluation of `f` with a (roughly) 
constant overhead, regardless of how many parameters `f` depends on.

!!! note
    The first time you run any gradient-based function (e.g. [`grad`](@ref)) there will
    likely be a delay of a few seconds while Enzyme does some precomputation and caching.

!!! warning
    The first time you run [`grad`](@ref) with a `General` mixer (including 
    [`mixer_clique`](@ref) and [`mixer_ring`](@ref)) you will likely get the warning:

    ```julia-repl
    Warning: Using fallback BLAS replacements, performance may be degraded
    ```

    This is a known [issue](https://discourse.julialang.org/t/warning-linking-two-modules-of-different-target-triples-bcloader-start/93864/4)
    with Enzyme, and should hopefully be addressed soon. It does not indicate any incorrect
    results, just that things are a bit slower than they need to be until some kinks between
    BLAS and Enzyme are worked out.

```@index
Pages = ["angle_finding.md"]
```

```@docs
grad(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
grad!(G::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals; flip_sign=false)
find_local_minimum(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
find_local_maximum(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
```
!!! note
    While many optimization packages, e.g. 
    [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html), 
    [`Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/), only include 
    minimization, we explicitly include both minimization and maximization. This is because 
    the standard trick of flipping the sign of the objective function to switch between 
    maximization and minimization doesn't work quite as smoothly with QAOA, as adding a 
    minus sign in front of ``H_C`` messess up the phases. Specifically, the fact that
    ```juli_repl
    exp_value(angles, mixer, obj_vals) != -exp_value(angles, mixer, -obj_vals)
    ```
    makes it so that ``\{\beta_i, \gamma_i\}`` which minimize a given QAOA do *not* maximize 
    that same QAOA under the replacement ``H_C \to -H_C``. In order to avoid having to 
    clarify whether some angles are good for e.g. "maximizing MaxCut" or 
    "minimizing -MaxCut", we provide separate implementations for maximization and 
    minimization which avoid this sign confusion.

```@docs
find_angles_bh(p, mixer::Mixer, obj_vals, measure=obj_vals; kwargs...)
grover_th(p, obj_vals; max=true)
```