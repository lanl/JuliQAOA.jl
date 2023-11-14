# Basic Use

The core functionality of JuliQAOA is to take in a set of angles 
``\{\beta_i, \gamma_i\}``, a mixer ``H_M``, and a cost function ``H_C``, and return the 
statevector

```math
|\psi(\beta, \gamma)\rangle = e^{-i \beta_p H_M} e^{-i \gamma_p H_C} \dots e^{-i \beta_1 
H_M} e^{-i \gamma_1 H_C} |\psi_0\rangle.
```

Here is a simple example for a 6-qubit MaxCut problem:

```julia
using JuliQAOA, Graphs

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
```

The statevector can then be used to calculate other quantities of interest, e.g. the 
expectation value of ``H_C`` or ground state probability. See the [Examples](@ref) section 
for more examples. 

### Angles
For a ``p``-round QAOA, the angles are given as a vector of length ``2p``, with the 
``\beta_i`` values followed by ``\gamma_i`` values. For example,
```julia
angles = [0.2, 3.1, 0.9, 1.4]
```
corresponds to a 2-round QAOA with 
``\beta_1 = 0.2, \beta_2 = 3.1, \gamma_1 = 0.9, \gamma_2 = 1.4``. 

### Mixer
JuliQAOA includes several built-in mixer types, along with the ability to define your own 
mixer. The standard transverse field mixer on `n` qubits is defined as `mixer_x(n)`, and 
returns an object of type `JuliQAOA.Mixer{X}` that stores all of the information about the 
mixer. Mixers are discussed more on the [Mixers](@ref) section. 

### Cost Function
Cost functions in JuliQAOA are given in terms of a function which takes as input an array of
0's and 1's and returns a real number. For example, one could define a cost function which 
measures how many adjacent bits are of opposite value like this:

```julia
function bitflips(x)
    C = 0
    for i in 1:length(x)-1
        if x[i] != x[i+1]
            C += 1
        end
    end
    return C
end
```

In order to simulate a QAOA with this cost function, simply evaluate it over all feasible 
states (we include iterators [`states(n)`](@ref) for all computational basis states, as well as [`dicke_states(n,k)`](@ref) for all states with `k` ones). This is commonly referred to as a list of objective values, which we abbreviate to 
`obj_vals`. In the case of unconstrained optimization problems, this is as straightforward 
as

```julia
obj_vals = [bitflips(x) for x in states(n)]
```

See the [Cost Functions](@ref) section for more information.




