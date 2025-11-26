# 3D Low-Frequency Magnetic FEM Solver (MATLAB)

A high-performance, fully parallelized finite element solver for low-frequency electromagnetics.

## Features
* **Formulation**: $\mathbf{A}-V$ potential formulation with Edge/Node mixed elements.
* **Solvers**: Static (Linear/Nonlinear), Transient (Eddy Current), Frequency Domain (AC).
* **Sources**: Meshless coil support via Biot-Savart integration.
* **HPC**: Fully parallelized assembly and field calculation (parfor).
* **Materials**: Support for nonlinear B-H curves and laminations.

## Quick Start

1.  **Setup**:
    Run `startup_fem.m` to add paths.

2.  **Run Examples**:
    * `solve_linear_static`: Basic iron core magnetization.
    * `solve_nonlinear_static`: Saturation test.
    * `solve_transient_eddy`: Skin effect step response.
    * `solve_frequency_eddy`: AC impedance extraction.

## Directory Structure
* `src/`: Core source code.
    * `assembly/`: Matrix assembly kernels.
    * `solver_*/`: Linear and nonlinear solvers.
    * `mesh_io/`: Gmsh file reader.
* `examples/`: Test scripts.

## Requirements
* MATLAB R2023b or newer.
* Parallel Computing Toolbox.
* Symbolic Math Toolbox (for kernel generation).
* (Optional) MUMPS solver interface.