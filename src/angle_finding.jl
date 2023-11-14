using Optim
using Basinhopping
using Enzyme
using LineSearches

"""
    grad(angles, mixer, obj_vals)

Calculate the gradient of [`exp_value(angles, mixer, obj_vals)`](@ref).

Can be extended to incorporate non-standard initial state ``|\\psi_0\\rangle`` and
observables in the same way as [`exp_value`](@ref).

!!! warning
    `grad` does not currently work with `Grover` mixers. For unconstrained
    problems you can use the equivalent `mixer_x(n, 0:n)/2^n`, and for constrained problems
    you can use `N = binom(n,k); mixer_general(dicke_states(n,k), ones(N,N)/N)`. These will
    both five you correct results, but are slower than the `Grover` implementation.
"""
function grad(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    G = similar(angles)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    grad!(G, sv, angles, mixer, obj_vals, measure)
    return G
end

function grad(sv::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    G = similar(angles)
    grad!(G, sv, angles, mixer, obj_vals, measure)
    return G
end

"""
    grad!(G, angles, mixer, obj_vals; flip_sign=false)

Calculate the gradient of [`exp_value(angles, mixer, obj_vals)`](@ref), storing the result
in the vector `G`. 

The optional argument `flip_sign` adds an overall minus sign to `G`, which can be necessary
to switch between using `grad!` to identify angles which minimize or maximize `exp_value`.

Can be extended to incorporate non-standard initial state ``|\\psi_0\\rangle`` and
observables in the same way as [`exp_value`](@ref).
"""
function grad!(G::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals; flip_sign=false)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    grad!(G, sv, angles, mixer, obj_vals, measure; flip_sign=flip_sign)
end

function grad!(G::Vector, sv::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals; flip_sign=false)
    sv_copy = copy(sv)
    dsv = zeros(ComplexF64, mixer.N)
    G .= 0.0
    f(a,b) = (flip_sign ? -1 : 1)*exp_value!(a, b, mixer, obj_vals, measure)
    Enzyme.autodiff(Reverse, f, Duplicated(sv_copy, dsv), Duplicated(angles, G))
end

"""
    find_local_minimum(angles, mixer, obj_vals)

Find ``\\{\\beta_i,\\gamma_i\\}`` which minimize the ``H_C`` defined by `obj_vals`,
beginning at `angles` and doing local search with 
[`BFGS`](https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm).

Angles are returned in the ranges ``\\beta_i \\in [0,T_{H_M}],~\\gamma_i \\in [0,T_{H_C}]``,
where the operator periods ``T`` are calculated with [`get_operator_period`](@ref).

Can be extended to incorporate non-standard initial state ``|\\psi_0\\rangle`` and
observables in the same way as [`exp_value`](@ref).
"""
function find_local_minimum(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    return find_local_minimum(sv, angles, mixer, obj_vals, measure)
end

function find_local_minimum(sv::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    ret = optimize(x->exp_value(sv, x, mixer, obj_vals, measure), 
                   (G,x)->grad!(G, sv, x, mixer, obj_vals, measure), 
                   angles, 
                   BFGS(linesearch=LineSearches.BackTracking()))
    return clean_angles(Optim.minimizer(ret), mixer, obj_vals)
end


"""
    find_local_maximum(angles, mixer, obj_vals)

Find ``\\{\\beta_i,\\gamma_i\\}`` which maximize the ``H_C`` defined by `obj_vals`,
beginning at `angles` and doing local search with 
[`BFGS`](https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm).

Angles are returned in the ranges ``\\beta_i \\in [0,T_{H_M}],~\\gamma_i \\in [0,T_{H_C}]``,
where the operator periods ``T`` are calculated with [`get_operator_period`](@ref).

Can be extended to incorporate non-standard initial state ``|\\psi_0\\rangle`` and
observables in the same way as [`exp_value`](@ref).
"""
function find_local_maximum(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    return find_local_maximum(sv, angles, mixer, obj_vals, measure)
end

function find_local_maximum(sv::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    ret = optimize(x->-exp_value(x, mixer, obj_vals, measure), 
                   (G,x)->grad!(G, sv, x, mixer, obj_vals, measure; flip_sign=true), 
                   angles, 
                   BFGS(linesearch=LineSearches.BackTracking()))
    return clean_angles(Optim.minimizer(ret), mixer, obj_vals)
end


"""
    get_operator_period(eigvals; tol=1e-10)

Calculate the period of the unitary transformation ``e^{i \\beta \\mathcal{O}}``, where the 
operator ``\\mathcal{O}`` has eigenvalues `eigvals`.

If `eigvals` is a list of integers, returns `2π/gcd(eigvals).` 
[`round`](https://docs.julialang.org/en/v1/base/math/#Base.round-Tuple{Type,%20Any})s 
entries of `eigvals` to integers within `tol`. Returns `Inf` if any of the `eigvals` are 
non-integers.
"""
function get_operator_period(eigvals; tol=1e-10)
    int_eigvals = round.(Int, eigvals)
    if all(abs.(eigvals .- int_eigvals) .< tol)
        return 2π/gcd(int_eigvals)
    else
        return Inf
    end
end

"""
    clean_angles(angles, mixer, obj_vals)

Takes `angles` and [`mod`](https://docs.julialang.org/en/v1/base/math/#Base.mod)s them into
the ranges ``\\beta_i \\in [0,T_{H_M}],~\\gamma_i \\in [0,T_{H_C}]``, where the operator
periods ``T`` are calculated with [`get_operator_period`](@ref).
"""
function clean_angles(angles, mixer, obj_vals)
    p = Int(length(angles)/2)
    if mixer.period < Inf
        bs = [mod(b, mixer.period) for b in angles[1:p]]
    else
        bs = angles[1:p]
    end
    obj_vals_period = get_operator_period(obj_vals)
    if obj_vals_period < Inf
        gs = [mod(g, obj_vals_period) for g in angles[p+1:end]]
    else
        gs = angles[p+1:end]
    end
    return [bs..., gs...]
end

function guess_angles(old_angles)
    if length(old_angles) == 0
        return rand(2)*2π
    else
        x0 = deepcopy(old_angles[end])
        rounds = Int(length(x0)/2)
        insert!(x0, rounds, x0[rounds])
        push!(x0, x0[end])
        return x0
    end
end

"""
    find_angles_bh(p, mixer, obj_vals; max=true, niter=100, file=nothing)

Find good angles up for the QAOA defined by `mixer`, `obj_vals` up to `p` rounds.

Uses an iterative, round-by-round angle finding algorithm that combines angle extrapolation
and [basinhopping](https://github.com/gamatos/Basinhopping.jl) to find high quality angles 
up to `p` rounds. See [this](https://arxiv.org/pdf/2202.00648.pdf) paper, section "Angle & 
Threshold Finding" for more details.

Optional arguments:
- `max=false`: determines whether the goal is to minimize or maximize `exp_value`
- `niter=100`: determines the number of basinhopping iterations
- `file=nothing`: save the resulting angles and expectation values in a plain text `file`
"""
function find_angles_bh(p, mixer::Mixer, obj_vals, measure=obj_vals; kwargs...)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    return find_angles_bh(sv, p, mixer, obj_vals, measure; kwargs...)
end

function find_angles_bh(sv::Vector, p::Int, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals; max=true, niter=100, file=nothing)
    # add code to load results if file exists
    # add code to check correct file type
    # add code to start skipping p if approximation ratio very close to 1
    angles = Vector{Vector{Float64}}()
    exp_vals = Vector{Float64}()
    opt = x0 -> optimize(x->(max ? -1 : 1)*exp_value(sv, x, mixer, obj_vals, measure), 
                   (G,x)->grad!(G, sv, x, mixer, obj_vals, measure; flip_sign=max), 
                   x0, 
                   BFGS(linesearch=LineSearches.BackTracking()))

    p = 1
    while i <= p
        x0 = guess_angles(angles)
        ret = basinhopping(opt, x0, BasinhoppingParams(niter=niter))
        new_angles = clean_angles(Optim.minimizer(ret), mixer, obj_vals)
        new_exp_val = (max ? -1 : 1)*minimum(ret)
        push!(angles, new_angles)
        push!(exp_vals, new_exp_val)
        println("done $i rounds, exp. value = $(round(new_exp_val, digits=5))")
        i += 1
        if !isnothing(file)
            save_angle_finding_results(file, angles, exp_vals)
        end
    end
    return angles, exp_vals
end

function save_angle_finding_results(file, angles, scores)
    open(file, "w") do f
        write(f, "angles = $(angles)\n")
        write(f, "scores = $(scores)")
    end
end