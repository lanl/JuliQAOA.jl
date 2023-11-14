# JuliQAOA

A fast, flexible package for simulating the Quantum Alternating Operator Ansatz (QAOA).

## Documentation

Please see our full documentation [here](https://lanl.github.io/JuliQAOA/stable/).

## Installation

The latest stable release of JuliQAOA can be installed using the Julia package manager with

```julia
] add https://github.com/lanl/JuliQAOA
```

## Usage

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
expectation value of ``H_C`` or ground state probability.

## Contributing

Please report any issues, bugs, feature requests, suggestions for improvement, etc., via the
Github **[issue tracker](https://github.com/lanl/JuliQAOA/issues)**. 

The primary developer of this package is [John Golden](mailto:golden@lanl.gov) 
([@johngolden](https://github.com/johngolden)). 

## License

This software is provided under a BSD license with a "modifications must be indicated"
clause. See the `LICENSE` file for the full text. 

**LANL C Number: C22038**

## Alternatives

QAOA can be simulated in general-purpose quantum simulators, e.g. 
[Qiskit](https://qiskit.org/documentation/stable/0.40/tutorials/algorithms/05_qaoa.html) 
and [Pennylane](https://pennylane.ai/qml/demos/tutorial_qaoa_intro/), however they will be
significantly slower.

[QAOA.jl](https://github.com/FZJ-PGI-12/QAOA.jl) is a circuit-based QAOA simulator for
Julia. [QOKit](https://github.com/jpmorganchase/QOKit/tree/main) is a Python package which
uses many of the same basic ideas as JuliQAOA, in particular precomputation and caching of
the cost function terms. It is currently more geared towards running highly parallelized 
simulations on large computer clusters.
