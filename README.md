# Heisenberg–Kitaev Simulations

This repository contains C and C++ codes I developed during my PhD in theoretical condensed matter physics to study disordered quantum spin systems, specifically the Heisenberg–Kitaev–Gamma model in high magnetic fields.

The goal of this repository is to provide examples of nontrivial scientific computing code written from scratch, illustrating numerical modelling, algorithmic implementation, and high-performance computation used throughout my doctoral research.

---

## Contents

### **1. `hgk_t-matrix/` – T-matrix calculation (C)**  
A standalone C implementation of a T-matrix calculation for single-impurity scattering.  
This code:
- encodes analytical expressions for Green’s functions,
- constructs and inverts matrices,
- computes impurity-induced corrections to observables.

**Compile:**  
```bash
gcc -O3 hkg_t-matrix.c -o hkg_t-matrix
```

### **2. `hkg_disorder-average/` – Real-space simulation with impurity averaging (C++)**
A C++ code performing real-space simulations of a disordered spin model with quenched bond disorder.
It:
-generates random bond disorder,
-performs iterative energy minimisation,
-computes observables,
-averages results over many disorder realisations,
-optionally parallelises this loop using OpenMP.

**Dependencies**

- C++ compiler (e.g. `g++`)
- **Eigen** (header-only linear algebra library): https://eigen.tuxfamily.org
- OpenMP (optional but recommended)

Eigen is not included in the repository.  
You must download it and ensure the `Eigen/` directory is visible to the compiler.

**Compilation**

If Eigen is installed system-wide or located in `/path/to/eigen`:

```bash
g++ -O3 -fopenmp -I /path/to/eigen hkg_disorder-average.cpp -o hkg_disorder-average
```

If OpenMP is unavailable, compile without it.
