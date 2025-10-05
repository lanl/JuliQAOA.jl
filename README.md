# MPS-JuliQAOA

> [!IMPORTANT]
> **Code updates:** We are currently updating the MPS code in this repository.  
> During this period, APIs, features, and file structures may change without notice.


A fast, flexible package for simulating the Quantum Alternating Operator Ansatz (QAOA) using the Matrix Product State (MPS).

## JuliQAOA

If you are interested in just JuliQAOA you can find the main branch [here](https://github.com/lanl/JuliQAOA.jl).

## Documentation

Please see our full documentation [here](https://lanl.github.io/JuliQAOA.jl/dev/).

## Activating an Environment

If you are new to Julia it is best practice to create a new environment for each project to avoid precompilation errors associated with version issues or package conflicts, and keep dependencies separated. This can be done with 

```julia
julia> # hit the `]` button to enter the package manager
(@v1.10) pkg> activate .
Activating project at `your:\folder\directory\project_name`
(project_name) pkg> activate .
```

## Installation

After you have activated your environment, the latest stable release of MPS-JuliQAOA can be installed using the Julia package manager with

```julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/lanl/JuliQAOA.jl#mps")
```
or 
```julia
julia> # hit the `]` button to enter the package manager
(project_name) pkg> add https://github.com/lanl/JuliQAOA.jl#mps
```

## Usage

The core functionality of MPS-JuliQAOA is to take in a set of angles 
${\beta_i, \gamma_i}$, a mixer $H_M$, and a cost function $H_C$, and return the 
statevector

```math
|\psi(\beta, \gamma)\rangle = e^{-i \beta_p H_M} e^{-i \gamma_p H_C} \dots e^{-i \beta_1 
H_M} e^{-i \gamma_1 H_C} |\psi_0\rangle.
```

Here is a simple example for a 3-qubit ZInteraction QAOA Hamiltonian:

```julia
using JuliQAOA
using Graphs, SimpleWeightedGraphs


#Create your own ZInteractions
interactions = [
    ZInteractions([1,2]),
    ZInteractions([2,3], 2)
]
```

Here is another example for which the helper function is used to build the MaxCut Hamiltonian:

```julia
using JuliQAOA
using Graphs, SimpleWeightedGraphs

#Custom graph
g = SimpleWeightedGraph(3)


add_edge!(g,1,2,1.0)
add_edge!(g,2,3,2.0)

interactions = maxcut_graph_to_zinteractions(g, weighted=true)
```

QAOA Problem:

```julia
using JuliQAOA, Graphs

#MPS-JuliQAOA code example...
```

Angle finding example:

```julia
using JuliQAOA, Graphs

#MPS-JuliQAOA code example...
```

The statevector can then be used to calculate other quantities of interest, e.g. the 
expectation value of ``H_C`` or ground state probability.

## Contributing

Please report any issues, bugs, feature requests, suggestions for improvement, etc., via the
Github **[issue tracker](https://github.com/lanl/JuliQAOA.jl/issues)**. 

The primary developer of this package is John Golden ([email](mailto:golden@lanl.gov), [github](https://github.com/johngolden)). 

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

## Citation

If you find MPS-JuliQAOA helpful in your work, please cite

```bibtex
@misc{feeney2025mpsjuliqaoauserfriendlyscalablempsbased,
      title={MPS-JuliQAOA: User-friendly, Scalable MPS-based Simulation for Quantum Optimization}, 
      author={Sean Feeney and Reuben Tate and John Golden and Stephan Eidenbenz},
      year={2025},
      eprint={2508.05883},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2508.05883}, 
}
```
