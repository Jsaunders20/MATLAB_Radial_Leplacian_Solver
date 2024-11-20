# Radial Polar Laplacian Solver in MATLAB

## Overview
This repository contains MATLAB code to solve a differential equation in polar coordinates using a finite difference method. The code aims to solve the Laplace equation for a two-dimensional domain with a radial component. Specifically, it computes the solution for a Laplacian in a polar grid using a finite difference scheme, with boundary conditions applied to simulate physical phenomena.

## Files
- **`radial_new.m`**: Sets up the polar grid, applies boundary conditions, and constructs the sparse matrix representation of the differential operator. The file also includes visualization of the boundary conditions applied to the grid.
- **`radialsolve_new.m`**: Solves the Laplace equation using the manufactured solution approach and computes the solution for the defined grid. It also provides a comparison with an exact solution.

## Methods and Techniques Used
### Finite Difference Scheme in Polar Coordinates
The main differential equation being solved is the Laplace equation in polar coordinates, expressed as:


$\frac{\partial^2 u}{\partial r^2} + \frac{1}{r} \frac{\partial u}{\partial r} + \frac{1}{r^2} \frac{\partial^2 u}{\partial \theta^2} = 0$


The solution approach uses a centered finite difference scheme to approximate the partial derivatives with respect to radial distance \(r\) and angular direction \(\theta\). The code handles both the inner and outer boundaries as well as periodic conditions in the angular direction.

### Sparse Matrix Representation
The Laplace operator is constructed using a sparse matrix representation (`A = sparse(N, N)`). The sparse matrix efficiently stores and computes the discretized Laplacian operator, which is necessary for handling large grid sizes. The entries of the matrix are populated based on the finite difference approximations.

### Boundary Conditions
The boundary conditions are defined by setting specific rows in the matrix `A` to represent known values on the boundaries:

- **Inner Radius Condition**: Enforces boundary conditions on the inner radius based on a sinusoidal function.
- **Outer Radius Condition**: Sets boundary conditions for the outer edge of the grid.
- **Angular Periodicity**: Handles the continuity of the solution in the angular direction by linking the first and last angular nodes.

Boundary conditions are visualized using a surface plot to show where they are applied.

## Details of the Code
### `radial_new.m`
This script:
1. Sets up the parameters and polar grid (θ, r).
2. Computes the sparse matrix `A` representing the Laplace operator using a finite difference scheme.
3. Initializes and visualizes boundary conditions using a manufactured boundary value function.
4. Constructs the boundary conditions matrix `DBDY` to indicate which nodes are boundary nodes and modifies `A` accordingly.
5. Visualizes the boundary condition application using `surf(X, Y, DBDY)`, where the `pol2cart` function is used to convert from polar to Cartesian coordinates for visualization purposes.

### `radialsolve_new.m`
This script:
1. Sets up a **manufactured solution** for consistency checking, with $u = r^2 \cos^2(θ)$.
2. Defines **boundary data** and source term `F` to solve the differential equation.
3. Computes the solution to the system `A \ rhs`, where `rhs` represents the right-hand side vector of boundary values and source terms.
4. Reshapes the solution vector back to the 2D polar grid for visualization.
5. Plots the numerical solution alongside the exact solution to illustrate the error and verify correctness.

## Running the Code
1. **Setup**: Ensure that MATLAB is installed and available in your environment.
2. **Execute**:
   - Run `radial_new.m` to create the grid, apply boundary conditions, and visualize the setup.
   - Run `radialsolve_new.m` to solve the Laplace equation, visualize the solution, and compare it to the exact solution.

## Visualization
- **Boundary Conditions**: A surface plot (`surf`) is used to show where the boundary conditions are applied to the domain, providing insight into the problem setup.
- **Solution Surface Plot**: The final solution is displayed in Cartesian coordinates with a comparison to the exact manufactured solution, highlighting areas where errors may occur.

## Notes
- The sparse matrix implementation significantly reduces memory consumption and computational time, making it feasible to solve larger grid sizes.
- Element-wise operations are indicated by the `.` operator (e.g., `.^`), which is crucial for correct matrix and vector operations.
- The `reshape` function is extensively used to transform between vector and matrix representations of the solution.

## Future Improvements
- Extend the solver to handle non-homogeneous boundary conditions or source terms.
- Optimize the grid refinement to achieve more accurate results with minimal computational cost.
- Implement additional visualization tools to better understand convergence and error distribution.

## Acknowledgements
This project serves as an educational tool for understanding finite difference methods in polar coordinates and sparse matrix manipulation in MATLAB. It can be a foundation for more complex numerical PDE solvers.

