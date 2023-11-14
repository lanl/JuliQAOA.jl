using JuliQAOA
using Test 

using LinearAlgebra
using Combinatorics
using Graphs

PX = [0 1; 1 0]
PY = [0 -im; im 0]
I = [1 0; 0 1]
H = [1 1; 1 -1]/sqrt(2)

function sum_symmetric_krons(arr)
    unique_perms = unique(collect(permutations(arr)))
    return sum([kron(x...) for x in unique_perms])
end

function qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)
    p = round(Int, length(angles)/2)
    sv = copy(init_state)
    for i in 1:p
        sv = exp.(-1im*angles[i+p]*obj_vals) .* sv
        sv = exp(-1im*angles[i]*mixing_mat) * sv
    end
    return sv
end

function random_sv(N)
    init_state = rand([-1,1],N).*rand(N) .+ rand([-1,1],N).*rand(N).*im
    init_state = init_state/norm(init_state)
    return init_state
end

@testset "JuliTest" begin

    include("utils.jl")

    include("cost_funcs.jl")

    include("mixers.jl")

    include("eval.jl")

    include("grover.jl")

    include("angle_finding.jl")
end
