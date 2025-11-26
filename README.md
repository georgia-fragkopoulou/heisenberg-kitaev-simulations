# Heisenberg–Kitaev Simulations

This repository contains C and C++ codes I developed during my PhD in theoretical condensed matter physics to study disordered quantum spin systems, in particular the Heisenberg–Kitaev–Gamma model in high magnetic fields.

The focus of these codes is on:
- Implementing analytical and numerical methods for quantum spin models.
- Performing real-space simulations with quenched disorder.
- Using high-performance computing concepts (parallelisation, efficient linear algebra) to obtain physically relevant observables.

The repository currently includes two main components:

1. **T-matrix calculation (`hkg_t-matrix/`)**  
   A C implementation of a T-matrix approach to study impurity scattering in a spin system. This code translates analytical expressions for Green’s functions and scattering amplitudes into a numerical procedure.

2. **Real-space impurity-averaged simulation (`hkg_disorder-average/`)**  
   A C++ code for real-space simulations of a disordered spin model with quenched bond disorder. It performs energy minimisation in the presence of impurities, diagonalizes the spin-wave Hamiltonian, and averages observables over many disorder realisations, using parallelisation for efficiency.

These codes produced numerical data used in my PhD thesis and related work on disordered quantum magnets.

> If you are reviewing this repository as part of an application: the goal here is to showcase examples of nontrivial scientific computing code written from scratch, rather than polished software packages.

## Contents

- `hkg_t-matrix/` – T-matrix calculation for a single impurity scattering.
- `hkg_disorder-average/` – Real-space simulation with disorder averaging.

## Dependencies and Compilation

The exact dependencies depend on the subproject; see the README in each subdirectory for details. In both cases, a standard C/C++ compiler (e.g. `gcc` or `g++`) is sufficient, with optional OpenMP support for parallelisation in the real-space code.

## Related Work

These codes were developed in the context of my PhD research on disordered quantum magnets. A more detailed description of the models and results can be found in:

- https://tud.qucosa.de/en/landing-page/https%3A%2F%2Ftud.qucosa.de%2Fapi%2Fqucosa%253A96696%2Fmets
- https://arxiv.org/abs/2408.12656
