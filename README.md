# Numerical Simulation of Heat Conduction

<div align="center">

![KHU Logo](https://theseedwikifile.theseed.io/7c/7ca289c1e4e3afb210607818d9c1ac79c13ad965f1a2ff30598d491b43452d4b.webp)

Department of Astronomy and Space Science, Kyung Hee University

Course: Numerical Simulation ([Prof. Tetsuya Magara](http://solardynamicslab.khu.ac.kr/~magara/))

</div>

---

## Overview

This project implements numerical solutions for the one-dimensional heat conduction equation as part of the Numerical Simulation course at the Department of Astronomy and Space Science, Kyung Hee University.

## Problem Description

The project solves the one-dimensional heat conduction equation:

<div align="center">

$$ \frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2} $$

</div>

where:

|  Symbol  | Description         |
| :------: | :------------------ |
|   $T$    | Temperature         |
|   $t$    | Time                |
|   $x$    | Position            |
| $\alpha$ | Thermal diffusivity |

## Implementation

### Numerical Method

- **FTCS** (Forward Time Central Space) explicit scheme
- Stability condition:
  <div align="center">

  $s = \frac{\alpha \Delta t}{(\Delta x)^2} \leq 0.5$

  </div>

### Code Structure

| File           | Description                                      |
| :------------- | :----------------------------------------------- |
| `heat.f90`     | Numerical solution only                          |
| `heat2.f90`    | Both numerical and analytical solutions          |
| `plotter-1.py` | Visualization for numerical solution             |
| `plotter-2.py` | Comparison of numerical and analytical solutions |

### Input Files

Parameters in `diff.in` / `diff2.in`:

|    Parameter    | Description                  |
| :-------------: | :--------------------------- |
| $\textit{jmax}$ | Number of grid points        |
| $\textit{nmax}$ | Maximum number of time steps |
|    $\alpha$     | Thermal diffusivity          |
|  $\textit{s}$   | Stability parameter          |
| $\textit{tmax}$ | Maximum simulation time      |

### Output Files

- Binary data: `diff.da` / `diff2.da`
- Text output: `diff.out` / `diff2.out`
- Visualization: Generated plots in `figs/` directory

## Usage

1. **Compile Fortran codes:**

   ```bash
   gfortran heat.f90 -o heat
   gfortran heat2.f90 -o heat2
   ```

2. **Run simulations:**

   ```bash
   ./heat    # Numerical solution
   ./heat2   # Numerical + Analytical solutions
   ```

3. **Generate plots:**
   ```bash
   python plotter-1.py   # Plot numerical solution
   python plotter-2.py   # Compare solutions
   ```

## Requirements

| Category  | Requirements                                    |
| :-------- | :---------------------------------------------- |
| Compiler  | Fortran compiler (gfortran recommended)         |
| Python    | Python 3.x                                      |
| Libraries | - NumPy<br>- Matplotlib<br>- Seaborn<br>- SciPy |
