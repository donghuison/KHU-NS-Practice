# Numerical Simulation of Heat Conduction

<div align="center">

![KHU Logo](https://theseedwikifile.theseed.io/7c/7ca289c1e4e3afb210607818d9c1ac79c13ad965f1a2ff30598d491b43452d4b.webp)

Department of Astronomy and Space Science, Kyung Hee University

Course: Numerical Simulation ([Prof. Tetsuya Magara](http://solardynamicslab.khu.ac.kr/~magara/))

[![Fortran](https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white)](https://fortran-lang.org/)
[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

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

### Numerical Discretization

Using the FTCS (Forward Time Central Space) scheme:

1. Time derivative (Forward difference):
   $$\frac{\partial T}{\partial t} \approx \frac{T_j^{n+1} - T_j^n}{\Delta t}$$

2. Space derivative (Central difference):
   $$\frac{\partial^2 T}{\partial x^2} \approx \frac{T_{j+1}^n - 2T_j^n + T_{j-1}^n}{(\Delta x)^2}$$

3. Resulting discrete equation:
   $$T_j^{n+1} = T_j^n + s(T_{j+1}^n - 2T_j^n + T_{j-1}^n)$$
   where $s = \frac{\alpha \Delta t}{(\Delta x)^2}$ is the stability parameter

### Stability Analysis

Von Neumann stability analysis shows that:
- The scheme is conditionally stable
- Stability condition: $s \leq 0.5$
- Physical interpretation: The time step must be small enough relative to the spatial discretization

### Analytical Solution

The analytical solution is derived using separation of variables and Fourier series expansion.

#### Boundary Value Problem

1. **Governing Equation:**
   $$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$

2. **Initial Condition:** Triangular temperature distribution
   $$T(x,0) = \begin{cases} 
   2x, & 0 \leq x < 0.5 \\
   2(1-x), & 0.5 \leq x \leq 1
   \end{cases}$$

3. **Boundary Conditions:** Fixed temperature at boundaries
   $$T(0,t) = T(1,t) = 0 \quad \text{for all } t \geq 0$$

#### Solution Process

1. **Step 1: Separation of Variables**
   - Assume solution of the form: $T(x,t) = X(x)G(t)$
   - Substituting into heat equation:
     $$X(x)G'(t) = \alpha X''(x)G(t)$$
   - Separating variables:
     $$\frac{G'(t)}{G(t)} = \alpha\frac{X''(x)}{X(x)} = -\alpha\lambda^2$$
   - This yields two ODEs:
     $$X''(x) + \lambda^2X(x) = 0$$
     $$G'(t) + \alpha\lambda^2G(t) = 0$$

2. **Step 2: Solving the ODEs**
   - Spatial solution: $X(x) = \sin(\lambda_n x)$, where $\lambda_n = (2n-1)\pi$
   - Temporal solution: $G(t) = e^{-\alpha \lambda_n^2 t}$

3. **Step 3: Fourier Series Expansion**
   - General solution: $T(x,t) = \sum_{n=1}^{\infty} a_n \sin(\lambda_n x)e^{-\alpha \lambda_n^2 t}$
   - Initial condition: 
   $$
   T(x,0) = \left\{
   \begin{aligned}
   2x     & \quad \text{for } 0 \leq x < 0.5 \\
   2(1-x) & \quad \text{for } 0.5 \leq x \leq 1
   \end{aligned}
   \right.
   $$
   - Fourier coefficients: $a_n = \frac{8}{(2n-1)^2\pi^2}(-1)^{\frac{2n-1-1}{2}}$


#### Numerical Implementation

```fortran
integer, parameter :: M_MAX = 100  ! Truncation of infinite series

do m=1, M_MAX
    dn = 2*m - 1                   ! Odd numbers for eigenvalues
    lamn = dn*pi                   ! Eigenvalues
    an = 8*(-1)**((dn-1)/2) / (dn**2 * pi**2)  ! Fourier coefficients
    dte = an*sin(lamn*x)*exp(-alpha*lamn**2*t) + dte
end do
```

#### Convergence Properties

1. **Series Truncation:**
   - The infinite series is truncated at M_MAX = 100 terms
   - Higher modes (n > 100) contribute negligibly due to:
     - Rapid decay of coefficients ($\sim \frac{1}{n^2}$)
     - Faster exponential decay for higher modes ($\sim e^{-n^2t}$)

2. **Physical Interpretation:**
   - Each term represents a standing wave mode
   - Lower modes persist longer (smaller $\lambda_n^2$ in exponential)
   - Initial sharp features require more modes for accurate representation
   - Solution smooths over time as higher modes decay rapidly


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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.