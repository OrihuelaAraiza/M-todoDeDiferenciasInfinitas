
# Finite Difference Solver

This project implements a finite difference solver for second-order differential equations. The solver uses finite differences to approximate the solution to differential equations under given boundary conditions.

* Valeria Liv Larumbe Colmenares
* Rodrigo López Moreno
* Jose Pablo Orihuela Araiza
* Luis Cedillo Maldonado

## Features

- Solve second-order differential equations using finite differences.
- Supports Dirichlet (value) and Neumann (derivative) boundary conditions.
- User-friendly CLI menu with options to run pre-defined test cases or solve custom equations.
- Verbose mode with pretty table output of results.
- Input validation with informative error messages.
- Colored output for improved readability using `colorama`.

## Dependencies

- Python 3
- `prettytable` for displaying results in a tabular format.
- `colorama` for colored CLI output.

To install the required dependencies, run:

```sh
pip install prettytable colorama
```

## How to Use

Run the script using Python:

```sh
python finite_difference_solver.py
```

Upon running the script, you will be presented with a menu:

1. **Run Test Cases**: Run a suite of pre-defined test cases to verify the solver's functionality.
2. **Solve Custom Equation**: Enter parameters for a custom differential equation, including boundary conditions and the function values.
3. **Exit**: Exit the program.

### Example Usage

**Test Case**

Select option `1` to run a simple test case with the following conditions:

- Differential Equation: \( y''(x) = 0 \)
- Interval: [0, 1]
- Boundary Conditions: \( y(0) = 0 \), \( y(1) = 1 \)

**Custom Equation**

Select option `2` and follow the prompts:

1. Enter the start and end of the interval (`a` and `b`), e.g., `a = 0`, `b = 1`.
2. Enter the number of internal points (`n`), e.g., `n = 10`.
3. Enter the function values for \( f(x) \), which defines the differential equation \( y''(x) = f(x) \).
4. Specify the type of boundary condition (`value` or `derivative`) and provide the appropriate value.
5. Confirm the displayed equation and boundary conditions before proceeding.

### Example Differential Equation

Consider the differential equation:

y''(x) = x

with the boundary conditions:

- Left Boundary Condition: \( y(0) = 1 \) (Dirichlet)
- Right Boundary Condition: \( y'(1) = -1 \) (Neumann)

Enter the corresponding values when prompted.

## Project Structure

- `finite_difference_solver.py`: Main script containing the solver implementation and CLI.
- `BoundaryCondition` class: Represents the boundary conditions (value or derivative).
- `FiniteDifferenceSolver` class: Implements the finite difference method to solve second-order differential equations.
