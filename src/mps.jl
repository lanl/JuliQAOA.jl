using ITensorMPS, ITensors

using Graphs, Random, CUDA, SimpleWeightedGraphs

using JuliQAOA
using BenchmarkTools

#Define operators that may be of use
ITensors.op(::OpName"ZZZ",::SiteType"Qubit") =
[1   0   0  0   0  0  0   0;
0  -1   0  0   0  0  0   0;
0   0  -1  0   0  0  0   0;
0   0   0  1   0  0  0   0;
0   0   0  0  -1  0  0   0;
0   0   0  0   0  1  0   0;
0   0   0  0   0  0  1   0;
0   0   0  0   0  0  0  -1;]

ITensors.op(::OpName"ZZ",::SiteType"Qubit") =
[1 0 0 0
0 -1 0 0
0 0 -1 0
0 0 0 1]

ITensors.op(::OpName"IDT",::SiteType"Qubit") =
[1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1]


ITensors.op(::OpName"X",::SiteType"Qubit") =
[0 1
1 0]

ITensors.op(::OpName"Z",::SiteType"Qubit") =
[1 0
0 -1]


struct ZInteractions #PAULIPOLY
    qubits::Vector{Int}     # qubit indices involved in the interaction
    weight::Float64       # coupling weight
end

ZInteractions(qubits::Vector{Int}) = ZInteractions(qubits, 1.0)

struct QAOAProblem
    interactions::Vector{ZInteractions}
    nqubits::Int
    sites::Vector{Index}
    psi0::MPS
end


function QAOAProblem(interactions::Vector{ZInteractions};
                          nqubits::Int = maximum([maximum(i.qubits) for i in interactions]),
                          site_type::String = "Qubit",
                          init_state::String = "+")  # could also be "Z0", "rand", etc.

    sites = siteinds(site_type, nqubits)

    # Default + state MPS
    psi0 = MPS(sites, init_state)

    return QAOAProblem(interactions, nqubits, sites, psi0)
end


function generate_phase_gates(problem::QAOAProblem, gamma::Float64)::Vector{ITensor}
    gates = ITensor[]

    for interaction in problem.interactions
        qubits = interaction.qubits
        weight = interaction.weight
        opstr = join(fill("Z", length(qubits)))
        #println("opstr $(opstr)")
        op_tensor = op(opstr, (problem.sites[q] for q in qubits)...)
        #@show norm(op_tensor)
        hj = weight * op_tensor
        Gj = exp(-im * gamma * hj)
        #@show Gj
        push!(gates, Gj)
    end

    return gates
end

function generate_xmix_gates(problem::QAOAProblem, beta::Float64)::Vector{ITensor}
    println("beta: $(beta)")
    Gj = [exp(-im * beta * (op("X", problem.sites[j]))) for j in 1:problem.nqubits]
    #println(Gj)
    return Gj
end

# A layer of the circuit we want to optimize
function layer(problem::QAOAProblem, angles::Vector)::Vector{ITensor}
    h_c_layer = generate_phase_gates(problem,angles[1])
    h_m_layer = generate_xmix_gates(problem,angles[2])
    return [h_c_layer; h_m_layer]
end


# The variational circuit we want to optimize
function variational_circuit(problem::QAOAProblem, angles::AbstractVector)::Vector{ITensor}

    p_rounds = Int(length(angles)/2)
    betas = angles[1:p_rounds]
    gammas = angles[p_rounds+1:2*p_rounds]
    circuit = layer(problem, [gammas[1],betas[1]])

    for i in 2:(p_rounds)
        println("round: $(i)")
        circuit = [circuit; layer(problem, [gammas[i],betas[i]])]
    end

    return circuit
end

function z_hamiltonian_opsum(interactions::Vector{ZInteractions})
    os = OpSum()
    for term in interactions
        qubits = term.qubits
        w = term.weight
        zops = [( "Z", q) for q in qubits]
        t = (w, [op for (op, q) in zops for op in (op, q)]...)

        os += t
        
    end
    #println("new os: $(os)")
    return os
end

function z_hamiltonian_mpo(problem::QAOAProblem)
    os = z_hamiltonian_opsum(problem.interactions)
    return MPO(os, problem.sites)
end



"""
Helper functions
"""
function maxcut_graph_to_zinteractions(g; weighted=false)
    if weighted != true
        interactions = ZInteractions[]
        for e in edges(g)
            x = src(e)
            y = dst(e)
            push!(interactions, ZInteractions([x, y]))
        end

    else
        interactions = ZInteractions[]
        weights(g)
        for e in edges(g)
            weight = get_weight(g, e)
            x = src(e)
            y = dst(e)
            push!(interactions, ZInteractions([x, y], weight))
        end
    end
    return interactions
end



function maxcut_post_process(e_val::Float64; interactions=nothing, weighted=false)

    le = length(interactions)#sum of weights for weighted graph
    post_e_val = le/2 - e_val/2

    return post_e_val
end

function mis_post_process(e_val::Float64; interactions=nothing, weighted=false, λ=2)

    #vertices/2 - (λ *numedges)/4 +e_val 
    #

    return post_e_val
end
"""
These are the gate based hamiltonians that are applie to the Matrix Product State (MPS)
"""



"""
The mixing Hamiltonian for QAOA gate based

"""


function run_qaoa_mps(angles::AbstractVector, problem::QAOAProblem; cutoff=1e-6, maxdim=64)
    H_mpo = z_hamiltonian_mpo(problem)
    #@show H_mpo
    U = variational_circuit(problem, angles)
    ψ_f = apply(U, problem.psi0; cutoff=cutoff, maxdim=maxdim)::MPS
    e_val = inner(ψ_f', H_mpo, ψ_f; cutoff=1e-4)
    #e_val = inner(problem.psi0', H_mpo, problem.psi0; cutoff=1e-4)
    return real(e_val)
end



function juli_qaoa(N::Integer, angles, g)
    # Generate random graph
    #Random.seed!(10)
    #g = random_regular_graph(N,3)#erdos_renyi(N, 0.5)
    obj_vals = [maxcut(g,x) for x in states(N)]
    mixer = mixer_x(N)
    x = exp_value(angles,mixer,obj_vals)
    return x
end



ITensors.op(::OpName"X",::SiteType"Qubit") =
[0 1
1 0]

ITensors.op(::OpName"Z",::SiteType"Qubit") =
[1 0
0 -1]

"""
GPU CODE

"""
function maxcut_ham_gpu(interactions::AbstractVector{<:Tuple{Vararg{Int}}}, sites::AbstractVector, gamma::Float64)::Vector{ITensor} 
    gates = ITensor[]
    #TODO @code_warntype macro says hj and Gj of type any, the inferring could hinder compilation performance
    for (i, j) in interactions
        hj = (0.5 * (op("IDT", sites[i],sites[j]) - op("ZZ", sites[i], sites[j])))
        hj = NDTensors.cpu(hj)
        Gj = exp(-im * gamma * hj)
        Gj = cu(Gj)
        push!(gates, Gj)
    end

    return gates
end

function x_mix_ham_gpu(N::Int64, sites::AbstractVector, beta::Float64)::Vector{ITensor}

    gates = ITensor[]
    
    for j in 1:N
        hj = op("X", sites[j])
        #exp() function performs operations which are not compatible with CUDA so take hj off GPU for exponential
        hj = NDTensors.cpu(hj)
        Gj = exp(-im * beta * hj)
        #load back onto GPU
        Gj = cu(Gj)
        push!(gates, Gj)
    end

    return gates
end
# The Matrix Product Operators of the Hamiltonians we are minimizing
function maxcut_ham_mpo(interactions)
    Z = [1 0 
         0 -1]

    ID = [1 0
          0 1]

    os = OpSum()

    for (i, j) in interactions

        os += 0.5, ID, i, ID, j
        os -= 0.5, "Z",  i, "Z",  j


    end
    #println("old os $(os)")
    return os
end

function spinglass_ham_mpo(instance::AbstractVector)#{<:Dict})
    
    os = OpSum()

    #linear terms
    for (key, val) in instance[1]
        os += val, "Z", key[1]
    end

    #quadratic terms
    for (key, val) in instance[2]
        os += val, "Z", key[1], "Z", key[2]
    end

    #quadratic terms
    for (key, val) in instance[3]
        os += val, "Z", key[1], "Z", key[2], "Z", key[3]
    end
    return os
end

"""
These are the gate based hamiltonians that are applie to the Matrix Product State (MPS)
"""

function maxcut_ham(interactions::AbstractVector, sites::AbstractVector, gamma::Float64)::Vector{ITensor} 

    Gj = [exp(-im * gamma * (0.5 * (op("IDT", sites[i],sites[j]) - op("ZZ", sites[i], sites[j])))) for (i, j) in interactions]
    println("OLD--------------------------")
    #@show Gj
    return Gj
end

function spinglass_ham(interactions::AbstractVector, sites::AbstractVector, gamma::Float64)::Vector{ITensor} 

    gates = ITensor[]

    #linear terms
    for (key, val) in interactions[1]
        hj = val * (op("Z", sites[key[1]]))
        Gj = exp(-im * gamma * hj)
        push!(gates, Gj)
    end

    #quadratic terms
    for (key, val) in interactions[2]
        hj = val * (op("ZZ", sites[key[1]], sites[key[2]]))
        Gj = exp(-im * gamma * hj)
        push!(gates, Gj)
    end

    #cubic terms
    for (key, val) in interactions[3]
        hj = val * (op( "ZZZ", sites[key[1]], sites[key[2]], sites[key[3]]))
        Gj = exp(-im * gamma * hj)
        push!(gates, Gj)

    end

    return gates
end


"""
The mixing Hamiltonian for QAOA gate based

"""

function x_mix_ham(N::Int64, sites::AbstractVector, beta::Float64)::Vector{ITensor}
    println("old beta: $(beta)")
    Gj = [exp(-im * beta * (op("X", sites[j]))) for j in 1:N]
    #println("Old Gj $(Gj)")
    return Gj
end

# A layer of the circuit we want to optimize
function layer(n::Int64, angles::Vector, sites::AbstractVector, interactions::AbstractVector)::Vector{ITensor}
    h_c_layer = maxcut_ham(interactions,sites,angles[1])
    h_m_layer = x_mix_ham(n,sites,angles[2])
    return [h_c_layer; h_m_layer]
    #return [h_m_layer;]

end
# A layer of the circuit we want to optimize
function layer_gpu(n::Int64, angles::Vector, sites::AbstractVector, interactions::AbstractVector)::Vector{ITensor}
    h_c_layer = maxcut_ham_gpu(interactions,sites,angles[1])
    h_m_layer = x_mix_ham_gpu(n,sites,angles[2])
    return [h_c_layer; h_m_layer]
end
# The variational circuit we want to optimize
function variational_circuit(n::Int64, angles::AbstractVector, sites::AbstractVector, interactions::AbstractVector)::Vector{ITensor}

    p_rounds = Int(length(angles)/2)
    betas = angles[1:p_rounds]
    gammas = angles[p_rounds+1:2*p_rounds]
    circuit = layer(n, [gammas[1],betas[1]], sites, interactions)

    for i in 2:(p_rounds)
        #println("round: $(i)")
        circuit = [circuit; layer(n, [gammas[i],betas[i]], sites, interactions)]
    end

    return circuit
 end

  # The variational circuit we want to optimize
function variational_circuit_gpu(n::Int64, angles::AbstractVector, sites::AbstractVector, interactions::AbstractVector)::Vector{ITensor}

    p_rounds = Int(length(angles)/2)
    betas = angles[1:p_rounds]
    gammas = angles[p_rounds+1:2*p_rounds]
    circuit = layer(n, [gammas[1],betas[1]], sites, interactions)

    for i in 2:(p_rounds)
        #println("round: $(i)")
        circuit = [circuit; layer(n, [gammas[i],betas[i]], sites, interactions)]
    end

    return circuit
end

function loss(angles::AbstractVector; kwargs...)
    H_opsum = maxcut_ham_mpo(kwargs[:interactions])
    H_mpo = MPO(H_opsum, kwargs[:sites])
    #@show H_mpo
    U = variational_circuit(kwargs[:N], angles, kwargs[:sites], kwargs[:interactions])
    ψ_f = apply(U, kwargs[:psi_0];  maxdim=kwargs[:bond_dimension], cutoff=1e-6)::MPS
    e_val = inner(ψ_f', H_mpo, ψ_f; cutoff=1e-4)
    #e_val = inner(kwargs[:psi_0]', H_mpo, kwargs[:psi_0]; cutoff=1e-4)
    return real(e_val)
end



function loss_gpu(angles::AbstractVector; kwargs...)
    H_opsum = maxcut_ham_mpo(kwargs[:interactions])
    H_mpo = MPO(H_opsum, kwargs[:sites])
    Cu_H = cu(H_mpo)
    U = variational_circuit(kwargs[:N], angles, kwargs[:sites], kwargs[:interactions])
    U = cu.(U)
    Cu_psi_0 = cu(kwargs[:psi_0])
    ψ_f = apply(U, Cu_psi_0; maxdim=kwargs[:bond_dimension])::MPS # cant use cutoff with GPU 
    x = inner(ψ_f', Cu_H, ψ_f)
    return -real(x)
end

function loss_cpu(angles::AbstractVector; kwargs...)
    H_opsum = maxcut_ham_mpo(kwargs[:interactions])
    H_mpo = MPO(H_opsum, kwargs[:sites])
    U = variational_circuit(kwargs[:N], angles, kwargs[:sites], kwargs[:interactions])
    ψ_f = apply(U, kwargs[:psi_0]; maxdim=kwargs[:bond_dimension], cutoff=1e-4)::MPS
    x = inner(ψ_f', H_mpo, ψ_f; cutoff=1e-4)
    return -real(x)
end

function loss_multithread_blas(angles::AbstractVector; kwargs...)
    H_opsum = maxcut_ham_mpo(kwargs[:interactions])
    H_mpo = MPO(H_opsum, kwargs[:sites])
    U = variational_circuit(kwargs[:N], angles, kwargs[:sites], kwargs[:interactions])
    ψ_f = apply(U, kwargs[:psi_0]; maxdim=kwargs[:bond_dimension], cutoff=1e-4)::MPS
    x = inner(ψ_f', H_mpo, ψ_f; cutoff=1e-4)
    return -real(x)
end

function loss_multithread_blocksparse(angles::AbstractVector; kwargs...)
    H_opsum = maxcut_ham_mpo(kwargs[:interactions])
    H_mpo = MPO(H_opsum, kwargs[:sites])
    U = variational_circuit(kwargs[:N], angles, kwargs[:sites], kwargs[:interactions])
    ψ_f = apply(U, kwargs[:psi_0]; maxdim=kwargs[:bond_dimension], cutoff=1e-4)::MPS
    x = inner(ψ_f', H_mpo, ψ_f; cutoff=1e-4)
    return -real(x)
end

function juli_qaoa(N::Integer, angles, g)
    obj_vals = [maxcut(g,x) for x in states(N)]
    mixer = mixer_x(N)
    x = exp_value(angles,mixer,obj_vals)
    return x
end
