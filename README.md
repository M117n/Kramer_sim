# Colloids Simulation Project

## Description

This repository contains the Fortran subroutines developed during a project focused on simulating colloidal systems under shear stress. The primary goal of the project was to investigate the response of a crystalline structure of colloids when subjected to various forces and stresses.

## Contents

- **`src/`**: Contains the main Fortran subroutines:
  - `GdR.f90`: Subroutine for calculating the radial distribution function (RDF) of the particle system.
  - `VERLET.f90`: Subroutine for generating the Verlet list to optimize force calculations between particles.

- **`data/`**: Placeholder for data files used in the simulations, such as initial particle conditions.

- **`images/`**: Placeholder for images and graphs generated from the simulation results.

- **`docs/`**: Additional documentation and the original project report in PDF format.

## Requirements

- Fortran compiler (e.g., `gfortran`)
- Basic knowledge of running Fortran programs

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/user/colloids-simulation.git
    cd colloids-simulation
    ```

2. Compile the Fortran subroutines:
    ```bash
    gfortran -o verlet src/VERLET.f90
    gfortran -o gdr src/GdR.f90
    ```

## Execution

To execute the subroutines, you can use the compiled binaries:

1. **Verlet List**:
    ```bash
    ./verlet
    ```
   This program generates the Verlet list, which is used to optimize the calculation of forces between particles within a certain cutoff radius.

2. **Radial Distribution Function (GdR)**:
    ```bash
    ./gdr
    ```
   This program calculates the radial distribution function (RDF), which helps in analyzing the spatial distribution of particles in the system.

## Results

Simulation results can be saved and analyzed. You can visualize them using external tools (like `gnuplot` or `matplotlib`).

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
