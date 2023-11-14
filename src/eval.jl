using LinearAlgebra


"""
    statevector(angles, mixer, obj_vals)

Calculate the QAOA statevector ``|\\psi(\\beta, \\gamma)\\rangle = e^{-i \\beta_p H_M} 
e^{-i \\gamma_p H_C} \\dots e^{-i \\beta_1 H_M} e^{-i \\gamma_1 H_C} |\\psi_0\\rangle``, 
where 
- `angles` = ``[\\beta_1, \\ldots, \\beta_p,\\gamma_1,\\ldots,\\gamma_p]``,
- `mixer` = ``H_M``,
- `obj_vals` = the diagonal of ``H_C``,
and ``|\\psi_0\\rangle`` is the uniform superposition over `mixer.feasible_states`.

    statevector(sv, angles, mixer, obj_vals)

Calculate the QAOA statevector with the custom initial state ``|\\psi_0\\rangle`` = `sv`.

# Examples

```julia-repl
using JuliQAOA, Graphs, LinearAlgebra

n = 6

# 3 rounds with random angles
p = 3 
angles = rand(2*p)

# transverse field mixer
mixer = mixer_x(n) 

# calculate the MaxCut cost function over all basis states on a random G(n,p) graph
g = erdos_renyi(n, 0.5)
obj_vals = [maxcut(g, x) for x in states(n)]

# calculate the statevector (with |ψ0⟩ = uniform superposition over all states)
statevector(angles, mixer, obj_vals)

# calculate the statevector (with |ψ0⟩ = random initial state)
sv = rand(2^n) + rand(2^n)*im
sv = sv/norm(sv)
statevector(sv, angles, mixer, obj_vals)
```
"""
function statevector(angles, mixer, obj_vals)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    statevector!(sv, angles, mixer, obj_vals)
    return sv
end

function statevector(sv, angles, mixer, obj_vals)
    sv2 = copy(sv)
    statevector!(sv2, angles, mixer, obj_vals)
    return sv2
end


"""
    statevector!(sv, angles, mixer, obj_vals)

Calculate the QAOA with ``|\\psi_0\\rangle`` = `sv`, and store the result in `sv`.
"""
function statevector!(sv, angles, mixer::Mixer{X}, obj_vals)
    p = Int(length(angles)/2)
    for i in 1:p
        applyExp!(sv, angles[i+p], obj_vals)
        applyH!(sv)
        applyExp!(sv, angles[i], mixer.d)
        applyH!(sv)
    end
end

function statevector!(sv, angles, mixer::Mixer{General}, obj_vals)
    p = Int(length(angles)/2)
    tmp = similar(sv)
    for i in 1:p
        applyExp!(sv, angles[i+p], obj_vals)
        mul!(tmp, mixer.vinv, sv)
        applyExp!(tmp, angles[i], mixer.d)
        mul!(sv, mixer.v, tmp)
    end
end


"""
    probabilities(angles, mixer, obj_vals)

Calculate the probability of observing each state in `mixer.feasible_states` for the QAOA
defined by `(angles, mixer, obj_vals)`, where``|\\psi_0\\rangle`` is the uniform 
superposition over `mixer.feasible_states`.

Equivalent to `abs2.(statevector(angles, mixer, obj_vals))`.

    probabilities(sv, angles, mixer, obj_vals)

Calculate the probability of observing each state in `mixer.feasible_states` for the QAOA
defined by `(angles, mixer, obj_vals)` with ``|\\psi_0\\rangle`` = `sv`.

Equivalent to `abs2.(statevector(sv. angles, mixer, obj_vals))`.
"""
function probabilities(angles, mixer, obj_vals)
sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
probabilities!(sv, angles, mixer, obj_vals)
return real.(sv)
end

function probabilities(sv, angles, mixer, obj_vals)
sv2 = copy(sv)
probabilities!(sv2, angles, mixer, obj_vals)
return real.(sv2)
end

"""
    probabilities!(sv, angles, mixer, obj_vals)

Calculate the probability of observing each state in `mixer.feasible_states` for the QAOA
defined by `(angles, mixer, obj_vals)` with ``|\\psi_0\\rangle`` = `sv`, storing the 
resulting probabilities in `sv`.
"""
function probabilities!(sv, angles, mixer, obj_vals)
    statevector!(sv, angles, mixer, obj_vals)
    for i in eachindex(sv)
        sv[i] = abs2(sv[i])
    end
end


"""
    exp_value(angles, mixer, obj_vals)

Calculate ``\\langle\\psi(\\beta, \\gamma)| H_C |\\psi(\\beta, \\gamma)\\rangle``, where
``|\\psi(\\beta, \\gamma)\\rangle`` = [`statevector(angles, mixer, obj_vals)`](@ref).

Equivalent to `dot(obj_vals, probabilities(angles, mixer, obj_vals))`.

Can be modified to measure the expectation value of something other than ``H_C`` and/or 
start at a non-standard ``|\\psi_0\\rangle`` as follows:

    exp_value(angles, mixer, obj_vals, measure)

Calculate ``\\langle\\psi(\\beta, \\gamma)| H_{\\text{measure}} |\\psi(\\beta, \\gamma)
\\rangle``, where `measure` = diagonal terms of ``H_{\\text{measure}}``.

Equivalent to `dot(measure, probabilities(angles, mixer, obj_vals))`.

    exp_value(sv, angles, mixer, obj_vals)

Calculate ``\\langle\\psi(\\beta, \\gamma)| H_C |\\psi(\\beta, \\gamma)\\rangle`` with 
``|\\psi_0\\rangle`` = `sv`.

    exp_value(sv, angles, mixer, obj_vals, measure)

Calculate ``\\langle\\psi(\\beta, \\gamma)| H_{\\text{measure}} |\\psi(\\beta, \\gamma)
\\rangle`` with ``|\\psi_0\\rangle`` = `sv`.


# Examples

```julia-repl
using JuliQAOA, Graphs, LinearAlgebra

n = 6

# 3 rounds with random angles
p = 3 
angles = rand(2*p)

# transverse field mixer
mixer = mixer_x(n) 

# calculate the MaxCut cost function over all basis states on a random G(n,p) graph
g = erdos_renyi(n, 0.5)
obj_vals = [maxcut(g, x) for x in states(n)]

# the traditional expectation value
exp_value(angles, mixer, obj_vals)

# the probability of observing an optimal state
measure = obj_vals .== maximum(obj_vals)
exp_value(angles, mixer, obj_vals, measure)

# the probability of observing an optimal state starting from a random initial state
sv = rand(2^n) + rand(2^n)*im
sv = sv/norm(sv)
exp_value(sv, angles, mixer, obj_vals, measure)
```
"""
function exp_value(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    sv = ones(ComplexF64, mixer.N)/sqrt(mixer.N)
    return exp_value!(sv, angles, mixer, obj_vals, measure)
end

function exp_value(sv::Vector, angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
    sv2 = copy(sv)
    return exp_value!(sv2, angles, mixer, obj_vals, measure)
end

"""
    exp_value!(sv, angles, mixer, obj_vals, measure)

Calculate ``\\langle\\psi(\\beta, \\gamma)| H_{\\text{measure}} |\\psi(\\beta, \\gamma)
\\rangle`` with ``|\\psi_0\\rangle`` = `sv`, storing the probabilities of observing each
state in `sv`.
"""
function exp_value!(sv, angles, mixer, obj_vals, measure)
    probabilities!(sv, angles, mixer, obj_vals)
    for i in eachindex(sv)
        sv[i] *= measure[i]
    end
    return real(sum(sv))
end

"""
Utilities
"""

"""
Calculate ``e^{-i a d}``, where `a` is a number and `d` is a vector, and store the result in
the vector `v`.
"""
function applyExp!(v, a, d)
    @inbounds for i in eachindex(v) 
        v[i] *= exp(-1im*a*d[i])
    end
end


"""
Apply the Hadamard gate ``H^{\\otimes n} to all qubits in `sv`, modifying `sv` in-place.
"""
function applyH!(sv)
    n = Int(log2(length(sv)))
    @inbounds for i in 1:n
        applyH!(sv,i)
    end
end

"""
Apply the Hadamard gate ``H^{\\otimes n} to the `l`th qubit in `sv`, modifying `sv` 
in-place.

Code modified (with permision) from https://blog.rogerluo.dev/2020/03/31/yany/.
"""
function applyH!(sv, l)
    U11 = 1/sqrt(2); U12 = 1/sqrt(2);
    U21 = 1/sqrt(2); U22 = -1/sqrt(2);
    step_1 = 1 << (l - 1)
    step_2 = 1 << l

    @inbounds if step_1 == 1
        for j in 0:step_2:size(sv, 1)-step_1
            ST1 = U11 * sv[j + 1] + U12 * sv[j + 1 + step_1]
            ST2 = U21 * sv[j + 1] + U22 * sv[j + 1 + step_1]

            sv[j + 1] = ST1
            sv[j + 1 + step_1] = ST2
        end
    elseif step_1 == 2
        for j in 0:step_2:size(sv, 1)-step_1
            Base.Cartesian.@nexprs 2 i->begin
                ST1 = U11 * sv[j + i] + U12 * sv[j + i + step_1]
                ST2 = U21 * sv[j + i] + U22 * sv[j + i + step_1]
                sv[j + i] = ST1
                sv[j + i + step_1] = ST2    
            end
        end
    elseif step_1 == 4
        for j in 0:step_2:size(sv, 1)-step_1
            Base.Cartesian.@nexprs 4 i->begin
                ST1 = U11 * sv[j + i] + U12 * sv[j + i + step_1]
                ST2 = U21 * sv[j + i] + U22 * sv[j + i + step_1]
                sv[j + i] = ST1
                sv[j + i + step_1] = ST2    
            end
        end
    elseif step_1 == 8
        for j in 0:step_2:size(sv, 1)-step_1
            Base.Cartesian.@nexprs 8 i->begin
                ST1 = U11 * sv[j + i] + U12 * sv[j + i + step_1]
                ST2 = U21 * sv[j + i] + U22 * sv[j + i + step_1]
                sv[j + i] = ST1
                sv[j + i + step_1] = ST2    
            end
        end
    else
        for j in 0:step_2:size(sv, 1)-step_1
            for i in j:8:j+step_1-1
                Base.Cartesian.@nexprs 8 k->begin
                    ST1 = U11 * sv[i + k] + U12 * sv[i + step_1 + k]
                    ST2 = U21 * sv[i + k] + U22 * sv[i + step_1 + k]
                    sv[i + k] = ST1
                    sv[i + step_1 + k] = ST2
                end
            end
        end
    end
end