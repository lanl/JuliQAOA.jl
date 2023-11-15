# JuliQAOA

A fast, flexible package for simulating the Quantum Alternating Operator Ansatz (QAOA).

## Documentation

Please see our full documentation [here](https://lanl.github.io/JuliQAOA.jl/dev/).

## Installation

The latest stable release of JuliQAOA can be installed using the Julia package manager with

```julia
] add https://github.com/lanl/JuliQAOA.jl
```

## Usage

The core functionality of JuliQAOA is to take in a set of angles 
${\beta_i, \gamma_i}$, a mixer $H_M$, and a cost function $H_C$, and return the 
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
# angles[1:p] = betas, angles[p+1:end] = gammas
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

If you find JuliQAOA helpful in your work, please cite

```bibtex
@inproceedings{10.1145/3624062.3624220, 
    author = {Golden, John and Baertschi, Andreas and O'Malley, Dan and Pelofske, Elijah and Eidenbenz, Stephan}, 
    title = {JuliQAOA: Fast, Flexible QAOA Simulation},
    year = {2023}, 
    isbn = {9798400707858},
    publisher = {Association for Computing Machinery}, 
    address = {New York, NY, USA}, 
    url = {https://doi.org/10.1145/3624062.3624220}, 
    doi = {10.1145/3624062.3624220},
    booktitle = {Proceedings of the SC '23 Workshops of The International Conference on High Performance Computing, Network, Storage, and Analysis}, 
    pages = {1454–1459}, 
    numpages = {6}, 
    location = {Denver, CO, USA}, 
    series = {SC-W '23} }
```
