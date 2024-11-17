# Lid Driven Cavity Flow
![C](https://img.shields.io/badge/language-C-00599C.svg) ![CMake](https://img.shields.io/badge/build-CMake-brightgreen.svg)

This repository contains a C implementation of a Computational Fluid Dynamics (CFD) simulation benchmark problem: the lid-driven cavity flow.

The lid-driven cavity flow problem is a classic CFD benchmark used to test the accuracy and efficiency of numerical solvers. The problem consists of a square cavity where the top lid moves with a constant velocity, while the other walls remain stationary. This creates a vortex flow pattern within the cavity that can be analyzed using numerical methods.

The results obtained using this code have been validated against the benchmark data provided in this paper:

> U. Ghia, K. N. Ghia, and C. T. Shin, "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method", *Journal of Computational Physics*, vol. 48, pp. 387–411, 1982.

**Note**: For a detailed comparison, please refer to the original paper [here](https://doi.org/10.1016/0021-9991(82)90058-4).

The implementation of this solver was inspired by a tutorial video that provided a clear explanation of CFD basics. You can watch the video [here](https://youtu.be/NWGYtdYCx3U) for more details.

## Requirements

- C compiler (e.g., `gcc`, `clang`, etc.)
- [CMake](https://cmake.org/download/) (version 3.10 or higher)
- Basic understanding of fluid dynamics and numerical methods.
- [git](https://git-scm.com/downloads) (optional, for version control)

## Getting Started

### Compilation

To compile the code, navigate to the project directory and use the following commands:

```bash
mkdir build
cd build
cmake ..
make
```
The compiled executable will be generated inside the `build` folder.

### Running the Simulation
Run the compiled executable:
```bash
./cfd_solver
```

### Output
The simulation generates a VTK file, which can be visualized using software like ParaView. The output file will be named `output.vtk` in the `build` folder.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.