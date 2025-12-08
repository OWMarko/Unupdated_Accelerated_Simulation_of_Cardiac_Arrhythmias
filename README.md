# Accelerated Simulation of Cardiac Arrhythmias in Heterogeneous Environments : A Projective Integration Approach (In progress since December 2025)

![Status](https://img.shields.io/badge/Status-Work_In_Progress-yellow)
![Language](https://img.shields.io/badge/C++-17-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## 📌 Overview
This repository hosts the source code for my **First year of MSc Thesis in Engineering Mathematics** (Université Côte d'Azur).

The goal is to develop a high-performance solver for cardiac electrophysiology using **Projective Integration** (Equation-Free framework). This method aims to overcome the time-scale stiffness inherent in biological models to accelerate "Digital Twin" simulations.

## The Mathematical Model
We focus on the **Aliev-Panfilov** model, a reaction-diffusion system describing the propagation of excitation waves in cardiac tissue:

$$
\frac{\partial u}{\partial t} = \nabla \cdot (\mathbf{D} \nabla u) - k u (u - a)(u - 1) - uv
$$

$$
\frac{\partial v}{\partial t} = \epsilon(u, v) (-v - k u (u - a - 1))
$$

Where $u$ is the transmembrane potential (fast variable) and $v$ is the recovery variable (slow variable).

## Architecture & Tech Stack
The project will be build to be modular and performance-oriented:

* **Core Solver :** Modern C++ (C++17)
* **Linear Algebra :** [Eigen](https://eigen.tuxfamily.org/) library for vectorized operations.
* **Parallelism :** OpenMP for multi-core spatial diffusion.
* **Numerical Scheme :** Hybrid approach combining an explicit micro-solver (Euler/RK4) and a projective macro-stepper.
* **Visualization :** VTK file export for analysis in [Paraview](https://www.paraview.org/).
* **Reproducibility : ** Docker containerization currently in the implementation pipeline

## Project Roadmap
This is a semester-long research project. Current progress:

- [x] **Phase 1:** Mathematical analysis of the Aliev-Panfilov system & stability conditions.
- [ ] **Literature Review:** In-depth study of selected papers on Projective Integration stability (In Progress).
- [ ] **Phase 2:** Repository setup and CMake architecture.
- [ ] **Phase 3:** Implementation of the 2D Laplacian operator using Finite Differences.
- [ ] **Phase 4:** Implementation of the Baseline Solver (Explicit Euler).
- [ ] **Phase 5:** Implementation of the **Projective Integration** algorithm (Time-jumps).
- [ ] **Phase 6:** Benchmarking & Paraview Visualization.

## Building the Project (Coming Soon)
Requirements: `CMake >= 3.10`, `g++` (supporting C++17), `Eigen3`.

```bash
mkdir build && cd build
cmake ..
make
./cardiac_solver
