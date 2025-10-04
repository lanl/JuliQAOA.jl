using ITensorMPS, ITensors

using Graphs, SimpleWeightedGraphs

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

function mis_graph_to_zinteractions()
    
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



