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

This project implements both numerical and analytical solutions for the one-dimensional heat conduction equation, developed as part of the Numerical Simulation course at the Department of Astronomy and Space Science, Kyung Hee University.

## Mathematical Foundation

### 1. Heat Conduction Equation

The project solves the one-dimensional heat conduction equation:

```math
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}
```

where:

|  Symbol  | Description         |
| :------: | :------------------ |
|   $T$    | Temperature         |
|   $t$    | Time                |
|   $x$    | Position            |
| $\alpha$ | Thermal diffusivity |

### 2. Numerical Solution (FTCS Scheme)

We use the FTCS (Forward Time Central Space) scheme for the numerical solution.

1. Time derivative (forward difference):
   $$\frac{\partial T}{\partial t} \approx \frac{T_j^{n+1} - T_j^n}{\Delta t}$$

2. Spatial second derivative (central difference):
   $$\frac{\partial^2 T}{\partial x^2} \approx \frac{T_{j+1}^n - 2T_j^n + T_{j-1}^n}{(\Delta x)^2}$$

3. Resulting discrete equation:
   $$T_j^{n+1} = T_j^n + s(T_{j+1}^n - 2T_j^n + T_{j-1}^n), \quad s = \frac{\alpha \Delta t}{(\Delta x)^2}.$$

#### Stability Analysis

Von Neumann stability analysis gives:

- The FTCS scheme is conditionally stable
- Stability condition: $s \leq 0.5$
- This means the time step must be sufficiently small relative to the spatial grid size.

### 3. Analytical Solution

The analytical solution is derived using separation of variables and Fourier series expansion.

#### Initial and Boundary Conditions

1. **Governing Equation:**

$$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$

2. **Initial Condition (Triangular profile):**

   This creates a symmetric “tent” shape peaking at $x=0.5$ with $T=1$ at the peak and $T=0$ at the boundaries.

$$
T(x,0) = \begin{cases}
2x,     & 0 \leq x < 0.5 \\
2(1-x), & 0.5 \leq x \leq 1
\end{cases} 
$$

3. **Boundary Conditions (Dirichlet):**

   $$T(0,t) = T(1,t) = 0 \quad \text{for all } t \geq 0$$

#### Solution Process

1. **Separation of Variables**

   - Assume solution of the form:
     $$T(x,t) = X(x)G(t)$$

   - Substituting into heat equation:
     <!-- $$X(x)G'(t) = \alpha X''(x)G(t)$$ -->
     $$X(x)G’(t) = \alpha X’’(x)G(t) \implies \frac{G’(t)}{G(t)} = \alpha \frac{X’’(x)}{X(x)} = -\alpha \lambda^2.$$

   <!-- - Separating variables:
     $$\frac{G'(t)}{G(t)} = \alpha\frac{X''(x)}{X(x)} = -\alpha\lambda^2$$ -->

   - This yields two ODEs:
     <!-- $$X''(x) + \lambda^2X(x) = 0 \quad \text{and} \quad G'(t) + \alpha\lambda^2G(t) = 0$$ -->
     $$X’’(x) + \lambda^2 X(x)=0, \quad G’(t) + \alpha \lambda^2 G(t)=0.$$


<!-- 2. **Step 2: Solving the ODEs** -->
2. **Eigenfunctions and eigenvalues**

   With boundary conditions $X(0)=0$, $X(1)=0$, we have:

   $$X_n(x)=\sin(n\pi x), \quad \lambda_n = n\pi, \quad n=1,2,3,\ldots.$$

   Due to the symmetry of the initial condition, only odd modes ($n=1,3,5,\ldots$) contribute. 

   By letting $n=2m-1$, we focus on these odd modes:

   $$\lambda_m = (2m-1)\pi, \quad X_m(x)=\sin((2m-1)\pi x).$$

   <!-- - Spatial solution:
     For the given boundary conditions $T(0,t) = T(1,t) = 0$, the spatial part of the solution satisfies a Sturm-Liouville problem. The general eigenfunctions are:
     $$X_n(x) = \sin(n \pi x), \quad n = 1, 2, 3, \ldots$$
     These eigenfunctions form a complete orthonormal basis on the interval $[0, 1]$ for functions that vanish at the boundaries. Each eigenfunction corresponds to an eigenvalue $\lambda_n = n \pi$.
     Due to the symmetry of the chosen initial condition (a triangular profile centered at $x = 0.5$), only odd modes $n = 1, 3, 5, \ldots$ will have nonzero Fourier coefficients. To simplify notation, we re-index the odd modes using $m$:
     $$n = 2m - 1, \quad m = 1, 2, 3, \ldots$$
     Thus, the eigenvalues and eigenfunctions relevant to this problem are:
     $$\lambda_m = (2m-1)\pi, \quad X_m(x) = \sin((2m-1)\pi x)$$ -->

3. **Temporal solution**
   For each mode:

   $$G_m(t) = e^{-\alpha (2m-1)^2 \pi^2 t}.$$
     <!-- The corresponding temporal equation $G'(t) + \alpha \lambda_m^2 G(t) = 0$ yields:
     $$G_m(t) = e^{-\alpha (2m-1)^2 \pi^2 t}$$ -->

#### Fourier Series Expansion

   - General solution:
     <!-- $$T(x,t) = \sum_{n=1}^{\infty} a_n \sin(\lambda_n x)e^{-\alpha \lambda_n^2 t}$$ -->
     <!-- Once the eigenfunctions and eigenvalues are established, the general solution can be written as a superposition of all modes:
     $$T(x,t) = \sum_{n=1}^{\infty} a_n \sin(\lambda_n x)e^{-\alpha \lambda_n^2 t}$$ -->
     $$T(x,t) = \sum_{m=1}^{\infty} a_m \sin((2m-1)\pi x) e^{-\alpha (2m-1)^2 \pi^2 t}.$$

   - Determining the coefficients $a_m$:

      - We project the initial condition onto the eigenfunctions:
      
      $$a_m = 2 \int_0^1 T(x,0) \sin((2m-1)\pi x),dx.$$
      
      - Performing this integral with the given piecewise linear initial condition yields:

      $$a_m = \frac{8}{(2m-1)^2 \pi^2}(-1)^{m-1}.$$
   
      The coefficients decrease as $1/(2m-1)^2$, ensuring rapid convergence, and the alternating sign $(-1)^{m-1}$ reflects the shape of the initial distribution.

      <!-- To find the coefficients $a_n$, we use the initial condition $T(x,0)$ and project it onto each eigenfunction $\sin(\lambda_n x)$:
      $$a_n = 2 \int_0^1 T(x,0) \sin(n \pi x) dx.$$
      Here, the factor of $2$ comes from the orthonormality of the sine functions under the given boundary conditions. Note that $\lambda_n = n \pi$. -->

   <!-- - Odd modes only:

      $$a_{2m-1} = 2 \int_0^1 T(x,0) \sin((2m-1) \pi x) dx, \quad m = 1, 2, 3, \ldots$$
      All other $a_n$ are zero. -->

   <!-- - Analytical form of $a_n$:

      Performing the integration for the given piecewise linear initial condition:

      ```math
         T(x,0) = \left\{
         \begin{aligned}
         2x,     & \quad 0 \leq x < 0.5 \\
         2(1-x), & \quad 0.5 \leq x \leq 1
         \end{aligned}
         \right.
      ```
      results in a closed-form expression for the coefficients. The final form for the nonzero coefficients $(odd \ n)$ is:
      $$a_{2m-1} = \frac{8}{(2m-1)^2\pi^2}(-1)^{m-1}$$
      This shows that the coefficients decreass as $1/n^2$, ensuring rapid convergence of the series.

   - Putting it all together:
     $$T(x,t) = \sum_{m=1}^{\infty} \left[ \frac{8}{(2m-1)^2\pi^2}(-1)^{m-1} \right] \sin((2m-1)\pi x)e^{-\alpha (2m-1)^2 \pi^2 t}.$$ -->

#### Numerical Implementation

An example Fortran snippet for reconstructing the analytical solution:
```fortran
integer, parameter :: M_MAX = 100  ! Truncation of infinite series
...
dte = 0.0
do m=1, M_MAX
    dn = 2*m - 1
    lamn = dn*pi
    an = 8.0*(-1)**(m-1) / ( (dn**2) * (pi**2) )
    dte = dte + an*sin(lamn*x)*exp(-alpha*(lamn**2)*t)
end do
```

#### Convergence Properties

1. **Series Truncation:**

   - Truncating at M_MAX=100 terms is sufficient because higher modes have negligible contribution.
   - Higher modes decay faster due to both the $1/n^2$ factor in coefficients and the exponential factor $e^{-n^2 t}$.

2. **Physical Interpretation:**
   - Each term represents a standing wave mode
   - Lower modes dominate the long-term solution; higher modes quickly vanish.
   - Over time, the temperature distribution becomes smoother as higher modes decay.

### Code Structure

| File           | Description                                      |
| :------------- | :----------------------------------------------- |
| `heat.f90`     | Numerical solution (FTCS) only                   |
| `heat2.f90`    | Numerical + analytical solutions                 |
| `plotter-1.py` | Visualization of the numerical solution          |
| `plotter-2.py` | Comparison of numerical and analytical           |

### Input Files

Parameters in `diff.in` / `diff2.in`:

|    Parameter    | Description                  |
| :-------------: | :--------------------------- |
| $\textit{jmax}$ | Number of grid points        |
| $\textit{nmax}$ | Max time steps               |
|    $\alpha$     | Thermal diffusivity          |
|  $\textit{s}$   | Stability parameter          |
| $\textit{tmax}$ | Maximum simulation time      |

### Output Files

- Binary data: `diff.da` / `diff2.da`
- Text output: `diff.out` / `diff2.out`
- Plots: Generated in `figs/` directory

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
   python plotter-2.py   # Compare numerical vs analytical solutions
   ```

## Requirements
   
| Category  | Requirements                                    |
| :-------- | :---------------------------------------------- |
| Compiler  | Fortran compiler (gfortran recommended)         |
| Python    | Python 3.x                                      |
| Libraries | - NumPy<br>- Matplotlib<br>- Seaborn<br>- SciPy |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
