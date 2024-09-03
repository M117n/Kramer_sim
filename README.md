# Colloids Simulation

 Description

This repository contains a Fortran-based Brownian Dynamics simulation. The primary goal of the project is to simulate the behavior of colloidal particles under various forces and stresses using the NVT ensemble and the Ermack & McCammon algorithm. The simulation includes calculations of radial distribution functions and uses the Verlet list for optimizing particle interaction computations.

## Contents

- **`src/`**: Contains the main Fortran program and subroutines:
  - `BD_v07.f90`: The main Brownian Dynamics simulation program.
  - `GdR.f90`: Subroutine for calculating the radial distribution function (RDF) of the particle system (integrated within `BD_v07.f90`).
  - `VERLET.f90`: Subroutine for generating the Verlet list to optimize force calculations between particles (integrated within `BD_v07.f90`).
  - `CORRELATION.f90`: Subroutine for computing the Mean Squared Displacement (MSD) as a correlation function over time (included in `BD_v07.f90`).

- **`data/`**: Placeholder for data files used in the simulations, such as initial particle configurations.

- **`docs/`**: Additional documentation, including the original project report in PDF format.

- **`LICENSE`**: License file for the project (LGPL-V3).

## Requirements

- Fortran compiler (e.g., `gfortran`)
- Basic knowledge of running Fortran programs

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/user/colloids-simulation.git
    cd colloids-simulation
    ```

2. Compile the Fortran program:
    ```bash
    gfortran -o brownian_dynamics src/BD_v07.f90
    ```

## Execution

To execute the Brownian Dynamics simulation, run the compiled program:

```bash
./brownian_dynamics
```

## Key Features:

- Verlet List: The program uses the Verlet list for efficient computation of particle interactions.
- Radial Distribution Function (GdR): The program computes the radial distribution function to analyze the spatial distribution of particles.
- Correlation Functions: The program computes the Mean Squared Displacement (MSD) as a correlation function over time using the CORRELATION subroutine.

## Simulation Results
Simulation results can be saved and analyzed. You can visualize them using external tools (like gnuplot or matplotlib). The program generates several output files, including:

- MovieBD.xyz: Particle positions.
- EPOT.txt: Potential energy.
- MSD.dat: Mean Squared Displacement (MSD) data.
