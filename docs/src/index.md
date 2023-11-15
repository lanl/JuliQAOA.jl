# JuliQAOA

A fast, flexible package for simulating the Quantum Alternating Operator Ansatz (QAOA).

Highlights include:
- No need for circuit-level description of QAOA problems.
- Exact statevector simulation with minimal memory overhead.
- Write cost functions directly in Julia, rather than encoding as a Hamiltonian.
- Native support for both constrained and unconstrained combinatorial optimization problems.
- Includes several built-in mixer Hamiltonians, user-defined mixers are also supported.
- Robust and extensible methods for learning good angles.

## Installation

The latest stable release of JuliQAOA can be installed using the Julia package manager with

```julia
import Pkg
Pkg.add(url="https://github.com/lanl/JuliQAOA.jl")
```
or 
```julia
] add https://github.com/lanl/JuliQAOA.jl
```

Test that the package works by running

```julia
] test JuliQAOA
```

## What is QAOA?

The Quantum Alternating Operator Ansatz (QAOA) is a quantum algorithm designed for solving 
combinatorial optimization problems. It is commonly studied as a hybrid classical/quantum 
heuristic, which combines quantum state evolution with classical parameter optimization. 
QAOA is based on the principle of evolving a quantum state through a series of operations to 
encode the solution to the optimization problem into the state of a quantum system.

The algorithm uses two Hamiltonians: ``H_C``, which encodes the cost function of the 
optimization problem, and ``H_M``, the mixing Hamiltonian, which facilitates transitions 
between different states. The goal is to prepare a quantum state whose measurement will 
reveal an optimal or near-optimal solution to the problem at hand. This is done by applying 
alternating unitary operators corresponding to the Hamiltonians. This evolution is dictated 
by the following sequence for ``p`` layers:

```math
|\psi(\beta, \gamma)\rangle = e^{-i \beta_p H_M} e^{-i \gamma_p H_C} \dots e^{-i \beta_1 
H_M} e^{-i \gamma_1 H_C} |\psi_0\rangle
```

The initial state ``|\psi_0\rangle`` is usually taken as the ground state of ``H_M``, which 
in many cases is the uniform superposition of feasible solutions. ``\beta_i`` and 
``\gamma_i`` are the variational parameters, often referred to as angles, that are optimized
classically. The objective is to choose them such that the expectation value of ``H_C`` is 
either minimized or maximized, depending on the optimization problem at hand:

```math
(\min/\max)_{\beta,\gamma} \langle\psi(\beta, \gamma)| H_C |\psi(\beta, \gamma)\rangle
```

Physically, one can think of tuning the angles to promote destructive interference between 
states with poor objective value and constructive interference between states with good 
objective value. QAOA can also be viewed as a Trotterization of quantum annealing.

QAOA was first introduced as the Quantum Approximate Optimization Algorithm by 
[Farhi et al.](https://arxiv.org/abs/1411.4028), and the framework was later generalized to 
the Quantum Alternating Operator Ansatz by 
[Hadfield et al.](https://arxiv.org/abs/1709.03489).

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
