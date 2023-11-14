using Memoize
using StatsBase: countmap

"""
Precompute and cache distinct objective values and their degeneracies.
"""
@memoize function initiate_grover(obj_vals)
    spectrum = countmap(obj_vals)
    scores = sort(collect(keys(spectrum)))
    degens = [spectrum[s] for s in scores]
    loc = Dict(x => i for (i, x) in enumerate(scores))
    return scores, degens, loc
end

"""
Calculate the expectation value with Grover mixer.
"""
function exp_value(angles::Vector, mixer::Mixer{Grover}, obj_vals::AbstractVector)
    scores, degens, _ = initiate_grover(obj_vals)
    coeffs = sv_grover(angles, scores, degens)
    for i in eachindex(coeffs)
        coeffs[i] = abs2(coeffs[i])*scores[i]*degens[i]
    end
    return real(sum(coeffs))
end

"""
Calculate the statevector value with Grover mixer.
"""
function statevector!(sv, angles, mixer::Mixer{Grover}, obj_vals)
    # Note: Grover mixer always works from uniform superposition
    # and therefore ignores the initial value of `sv`
    scores, degens, loc = initiate_grover(obj_vals)
    coeffs = sv_grover(angles, scores, degens)
    for i in eachindex(sv)
        sv[i] = coeffs[loc[obj_vals[i]]]
    end
end

"""
Base functionality for the efficient Grover statevector simulation.
"""
function sv_grover(angles, scores, degens)
    N = sum(degens)
    coeffs = ones(ComplexF64, length(scores))/sqrt(N)
    p = Int(length(angles)/2)
    for i in 1:p
        applyExp!(coeffs, angles[i+p], scores)
        state_sum = sum(coeffs .* degens)*(1-exp(-1im*angles[i]))/N
        coeffs .-= state_sum
    end
    return coeffs
end

"""
    grover_th_ev(p, obj_vals; max=true)

Returns the optimal expectation value for the Grover-Th QAOA variant introduced
[here](https://arxiv.org/pdf/2106.13860.pdf). 

Grover-Th represent a direct port of Grover's search algorithm in the QAOA context. In other 
words, the output of this function represents the expectation value at `p` rounds that one
could recover by simply doing unstructured search for optimal states (technically it's a bit 
more complicated than this, see [here](https://arxiv.org/pdf/2202.00648.pdf) for a more 
precise characterization). This can serve as a nice benchmark for QAOA performance, that is,
QAOA should at least be able to beat unstructured search.

Set `max=false` for minimization problems.
"""
function grover_th(p, obj_vals; max=true)
    th = optimal_grover_th(p, obj_vals; max=max)
    return exp_value_grover_th(p, obj_vals, th; max=max)
end

"""
Find the optimal threshold to set at `p` rounds to maximize (or minimize) the expectation
value. 
"""
function optimal_grover_th(p, obj_vals; max=true)
    scores, _, _ = initiate_grover(obj_vals)
    best_ev = 0
    best_th = 0
    for th in scores
        ev = exp_value_grover_th(p, obj_vals, th; max=max)
        if max ? ev > best_ev : ev < best_ev
            best_ev = ev
            best_th = th
        end
    end
    return best_th
end

"""
Calculate the expectation value for Grover-Th.
"""
function exp_value_grover_th(p, obj_vals, th; max=true)
    N = length(obj_vals)
    scores, degens, _ = initiate_grover(obj_vals)
    num_fail_th = 0
    scores_th = similar(scores)
    for i in 1:length(scores)
        if max ? scores[i] <= th : scores[i] >= th
            num_fail_th += degens[i]
            scores_th[i] = 0
        else
            scores_th[i] = 1
        end
    end
    r = num_fail_th/N
    angles = optimal_grover_angles(p, N, r)
    coeffs = sv_grover(angles, scores_th, degens)
    for i in eachindex(coeffs)
        coeffs[i] = abs2(coeffs[i])*scores[i]*degens[i]
    end
    return real(sum(coeffs))
end

"""
Calculate the optimal angles for Grover-Th.
"""
function optimal_grover_angles(p, obj_vals::AbstractArray, th::Number; max=true)
    N = length(obj_vals)
    scores, degens, _ = initiate_grover(obj_vals)
    num_fail_th = 0
    scores_th = similar(scores)
    for i in 1:length(scores)
        if max ? scores[i] <= th : scores[i] >= th
            num_fail_th += degens[i]
            scores_th[i] = 0
        else
            scores_th[i] = 1
        end
    end
    r = num_fail_th/N
    return optimal_grover_angles(p, N, r)
end

function optimal_grover_angles(p, N::Number, r::Number)
    c0 = 1/sqrt(N) # coefficient of states passing the threshold
    c1 = 1/sqrt(N) # coefficient of states failing the threshold
    bs = Array{Float64,1}(); gs = Array{Float64,1}()
    while length(bs) < p && 4*(1-r)/N < c0^2
        state_sum = 2*(r*c0-(1-r)*c1)
        c0 = c0 - state_sum
        c1 = -c1 - state_sum
        push!(bs, pi)
        push!(gs, pi)
    end
    if length(bs) < p && 4*(1-r)/N >= c0^2
        if r == 1/2
            push!(bs,pi/2)
            push!(gs,pi/2)
        else
            beta_opt_x = -abs(c0)*sqrt(4*(1-r)/N-c0^2)
            beta_opt_y = 2*(1-r)/N - c0^2
            
            gamma_opt_x = -sqrt(4*(1-r)/N-c0^2)/(c1*sign(c0))
            gamma_opt_y = c0*(1-2*r)/c1
            
            beta_opt = atan(beta_opt_x,beta_opt_y)
            gamma_opt = atan(gamma_opt_x,gamma_opt_y)
            
            push!(bs, beta_opt)
            push!(gs, gamma_opt)
        end
    end
    while length(bs) < p
        push!(bs,0)
        push!(gs,0)
    end
    return [bs..., gs...]
end
