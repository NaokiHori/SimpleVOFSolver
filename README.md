# Simple VOF Solver

[![License](https://img.shields.io/github/license/NaokiHori/SimpleVOFSolver)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml)
[![DOC](https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/documentation.yml/badge.svg)](https://naokihori.github.io/SimpleVOFSolver/index.html)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SimpleVOFSolver/main)](https://github.com/NaokiHori/SimpleVOFSolver/commits/main)

[![Simulation Snapshot](https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot2d.png)](https://youtu.be/CpqQJxSkm7Q)

## Overview

This library numerically simulates the motion of free surfaces and deformable droplets in two- and three-dimensional Cartesian domains using the finite-difference and volume-of-fluid methods. It is built on top of [`SimpleNSSolver`](https://github.com/NaokiHori/SimpleNSSolver).

## Dependency

This solver shares the same dependencies as [`SimpleNSSolver`](https://github.com/NaokiHori/SimpleNSSolver).

## Quick Start

1. **Prepare the workspace**

   ```sh
   mkdir -p /path/to/your/directory
   cd /path/to/your/directory
   ```

2. **Clone the repository**

   ```sh
   git clone --recurse-submodules https://github.com/NaokiHori/SimpleVOFSolver
   cd SimpleVOFSolver
   ```

3. **Set initial conditions**

   `Python3` is used to initialize the flow fields conveniently. You can provide `NPY` files differently under `initial_condition/output/`.

   ```sh
   cd initial_condition
   make output
   bash exec.sh
   cd ..
   ```

4. **Build the solver**

   ```sh
   make output
   make all
   ```

5. **Run the simulation**

   ```sh
   bash exec.sh
   ```

## Documentation

A brief overview of the free-surface treatment is available [here](https://naokihori.github.io/SimpleVOFSolver). For additional details, please refer to the [`SimpleNSSolver` documentation](https://naokihori.github.io/SimpleNSSolver).

## Acknowledgements

The volume-of-fluid method implemented in this solver is based on the [`THINC/QQ scheme`](https://www.sciencedirect.com/science/article/pii/S0021999117305995?via%3Dihub) with some modifications.

