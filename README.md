# 💎 JEWEL.jl

> **Julia Elastic Wave Equation Laboratory**
>
> From equations to earthquakes — learn 3D seismic wave propagation one line of Julia at a time.

---

## Overview

JEWEL (Julia Elastic Wave Equation Laboratory) is an educational implementation of a three-dimensional Finite Difference Time Domain (FDTD) solver for the elastic wave equation, written entirely in Julia.

Unlike most wave propagation codes developed for large-scale HPC simulations, JEWEL focuses on readability, transparency, and learning. Every numerical operation is intentionally written in a clear and explicit way so that students and researchers can understand how a wave equation solver is built from scratch.

The project aims to bridge the gap between theoretical seismology textbooks and production-scale numerical simulators.

---

## Philosophy

Most numerical wave propagation software is designed to answer research questions.

JEWEL is designed to answer questions like:

• How is the elastic wave equation discretized?
• How are finite differences implemented?
• How are absorbing boundary conditions applied?
• How does a moment tensor source generate seismic waves?
• How are displacement, velocity, and acceleration computed?
• How are synthetic seismograms recorded?

If you've ever wanted to understand every line of a seismic wave propagation code, this project is for you.

---

## Features

✓ 3D elastic wave propagation

✓ Explicit Finite Difference Time Domain (FDTD)

✓ Isotropic elastic medium

✓ Constant material properties

✓ Force source implementation

✓ Double-couple moment tensor source

✓ Gaussian and Ricker source wavelets

✓ Absorbing and reflecting boundary conditions

✓ Receiver recording (displacement, velocity and acceleration)

✓ Wavefield visualization using Makie

✓ Automatic MP4 animation generation

✓ Clean and well-commented Julia implementation

---

## Numerical Method

Current implementation:

• Explicit second-order time integration
• Centered finite-difference spatial derivatives
• Isotropic elastic wave equation
• Displacement formulation
• Cerjan-style absorbing boundaries
• Moment tensor and force sources

Future versions will include

• Variable velocity models
• Free surface boundary conditions
• Perfectly Matched Layers (PML)
• Higher-order finite differences
• Anisotropy
• GPU acceleration
• Parallel computing
• Topography
• Viscoelastic attenuation

---

## Repository Structure

JEWEL/

├── src/

│   ├── JEWEL_functions.jl

│   ├── JEWEL_structures.jl

│   └── derivative_constants.jl

│

├── pics/

│

├── RUN_JEWEL.jl

│

└── README.md

---

## Getting Started

Clone the repository

git clone https://github.com/yourusername/JEWEL.jl

Enter the project

cd JEWEL.jl

Install the required packages

using Pkg
Pkg.instantiate()

Run the example

include("RUN_JEWEL.jl")

---

## Example Workflow

1. Define model dimensions.

2. Define elastic properties (Vp, Vs and density).

3. Choose a source type:
   - Point force
   - Moment tensor

4. Define boundary conditions.

5. Place receivers.

6. Run the simulation.

7. Visualize the propagating wavefield.

8. Export synthetic seismograms.

9. Create an MP4 animation.

---

## Educational Roadmap

Part I
Constant isotropic elastic medium

Part II
Finite-difference derivatives

Part III
Boundary conditions

Part IV
Seismic source implementation

Part V
Receiver recording

Part VI
Attenuation

Part VII
Variable elastic models

Part VIII
Free surface

Part IX
GPU implementation

Each chapter introduces one new concept while keeping the numerical implementation easy to follow.

---

## Why JEWEL?

Modern wave propagation software can contain hundreds of thousands of lines of code.

JEWEL intentionally avoids unnecessary complexity.

Instead of asking

"How fast can we simulate?"

it asks

"Can someone understand every line?"

The goal is not to compete with production software.

The goal is to help build the next generation of computational seismologists.

---

## Planned Features

□ Heterogeneous models

□ PML absorbing boundaries

□ Free surface

□ Anisotropic elasticity

□ Multiple simultaneous sources

□ DAS (Distributed Acoustic Sensing) receivers

□ GPU support

□ MPI parallelization

□ Adaptive visualization

□ Interactive notebooks

---

## Citation

If JEWEL contributes to your research or teaching, please consider citing the repository once a DOI becomes available.

---

## Contributing

Contributions are welcome.

Whether it is

• fixing bugs

• improving documentation

• adding numerical methods

• implementing new boundary conditions

• optimizing performance

or simply improving readability,

every contribution helps make JEWEL a better educational resource.

---

## License

MIT License

---

## Author

Developed by Mariano Simón Arnaiz Rodríguez

Institut de Physique du Globe de Paris (IPGP)

Julia Elastic Wave Equation Laboratory

Making computational seismology understandable.
