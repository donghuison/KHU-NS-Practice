# Numerical Simulation of Heat Conduction

## Overview

This project implements numerical solutions for the one-dimensional heat conduction equation as part of the Numerical Simulation course at the Department of Astronomy and Space Science, Kyung Hee University.

## Problem Description

The project solves the one-dimensional heat conduction equation:

$$ \frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2} $$

where:

- $T$: Temperature
- $t$: Time
- $x$: Position
- $\alpha$: Thermal diffusivity

## Implementation

### Numerical Method

- FTCS (Forward Time Central Space) explicit scheme
- Stability condition: $s = αΔt/Δx² ≤ 0.5$

### Code Structure

1. `heat.f90`: Numerical solution only
2. `heat2.f90`: Both numerical and analytical solutions
3. `plotter-1.py`: Visualization for numerical solution
4. `plotter-2.py`: Comparison of numerical and analytical solutions

### Input Files

- `diff.in` / `diff2.in`: Contains simulation parameters
  ```
  jmax nmax alpha s tmax
  ```
  - jmax: Number of grid points
  - nmax: Maximum number of time steps
  - alpha: Thermal diffusivity
  - s: Stability parameter
  - tmax: Maximum simulation time

### Output Files

- `diff.da` / `diff2.da`: Binary data files
- `diff.out` / `diff2.out`: Text output files
- Generated plots in `figs/` directory

## Usage

1. Compile Fortran codes:
   ```bash
   gfortran heat.f90 -o heat
   gfortran heat2.f90 -o heat2
   ```
2. Run simulations:
   ```bash
   ./heat    # Numerical solution
   ./heat2   # Numerical + Analytical solutions
   ```
3. Generate plots:
   ```bash
   python plotter-1.py   # Plot numerical solution
   python plotter-2.py   # Compare solutions
   ```

## Requirements

- Fortran compiler (gfortran compiler recommended)
- Python 3.x
- NumPy
- Matplotlib
- Seaborn
- SciPy
