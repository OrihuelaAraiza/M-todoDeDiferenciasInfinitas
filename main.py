import time
from fractions import Fraction
from prettytable import PrettyTable
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)


class BoundaryCondition:
    """
    Represents a boundary condition for the solver.
    Can be a value (Dirichlet) or a derivative (Neumann).
    """

    def __init__(self, value=None, derivative=None):
        if value is not None and derivative is not None:
            raise ValueError("A boundary condition cannot have both value and derivative.")
        self.value = value
        self.derivative = derivative


class FiniteDifferenceSolver:
    """
    Solver for second-order differential equations using finite differences.
    """

    EPSILON = 1e-12  # Tolerance to avoid division by zero

    def __init__(self, func, a, b, n, left_bc, right_bc):
        """
        Initializes the solver.
        :param func: Callable function f(x).
        :param a: Start of the interval (float).
        :param b: End of the interval (float).
        :param n: Number of internal points (int).
        :param left_bc: BoundaryCondition at x=a.
        :param right_bc: BoundaryCondition at x=b.
        """
        if not callable(func):
            raise ValueError("The function 'func' must be callable.")
        if a >= b:
            raise ValueError("The interval start 'a' must be less than 'b'.")
        if n < 2:
            raise ValueError("The number of points 'n' must be at least 2.")
        self.func = func
        self.a = a
        self.b = b
        self.n = n
        self.left_bc = left_bc
        self.right_bc = right_bc
        self.h = (b - a) / n
        self.x_points = [a + i * self.h for i in range(n + 1)]

    def _create_matrix_and_vector(self):
        """
        Creates the diagonals for the tridiagonal system and the vector b.
        :return: main_diag, lower_diag, upper_diag, b
        """
        size = self.n + 1
        main_diag = [-2 for _ in range(size)]
        lower_diag = [1 for _ in range(size - 1)]
        upper_diag = [1 for _ in range(size - 1)]
        b = [self.h**2 * self.func(self.x_points[i]) for i in range(size)]

        # Apply boundary conditions
        self._apply_boundary_conditions(main_diag, lower_diag, upper_diag, b)

        return main_diag, lower_diag, upper_diag, b

    def _apply_boundary_conditions(self, main_diag, lower_diag, upper_diag, b):
        """
        Applies the boundary conditions to the diagonals and vector b.
        """
        # Left boundary condition
        if self.left_bc.value is not None:  # Dirichlet
            main_diag[0] = 1
            b[0] = self.left_bc.value
            upper_diag[0] = 0
        elif self.left_bc.derivative is not None:  # Neumann
            main_diag[0] = -1 / self.h
            upper_diag[0] = 1 / self.h
            b[0] = self.left_bc.derivative

        # Right boundary condition
        if self.right_bc.value is not None:  # Dirichlet
            main_diag[-1] = 1
            b[-1] = self.right_bc.value
            lower_diag[-1] = 0
        elif self.right_bc.derivative is not None:  # Neumann
            lower_diag[-1] = -1 / self.h
            main_diag[-1] = 1 / self.h
            b[-1] = self.right_bc.derivative

    def _solve_tridiagonal_system(self, main_diag, lower_diag, upper_diag, b):
        """
        Solves the tridiagonal system using forward and back substitution.
        :param main_diag: Main diagonal of the matrix.
        :param lower_diag: Lower diagonal of the matrix.
        :param upper_diag: Upper diagonal of the matrix.
        :param b: Right-hand side vector.
        :return: Solution vector y.
        """
        n = len(main_diag)

        # Forward sweep
        for i in range(1, n):
            factor = lower_diag[i - 1] / main_diag[i - 1]
            main_diag[i] -= factor * upper_diag[i - 1]
            b[i] -= factor * b[i - 1]

        # Back substitution
        y = [0 for _ in range(n)]
        y[-1] = b[-1] / main_diag[-1]
        for i in range(n - 2, -1, -1):
            y[i] = (b[i] - upper_diag[i] * y[i + 1]) / main_diag[i]

        return y

    def solve(self, verbose=False):
        """
        Solves the differential equation using finite differences.
        :param verbose: If True, prints step-by-step solution.
        :return: List of y values at each x point.
        """
        start_time = time.time()
        main_diag, lower_diag, upper_diag, b = self._create_matrix_and_vector()
        solution = self._solve_tridiagonal_system(main_diag, lower_diag, upper_diag, b)
        end_time = time.time()

        if verbose:
            print(Fore.CYAN + f"\nSolver completed in {end_time - start_time:.4f} seconds.")
            table = PrettyTable()
            table.field_names = ["x", "y"]
            for i, y_val in enumerate(solution):
                x_str = f"{self.x_points[i]:.2f}"
                if abs(y_val) < 1e-6:
                    y_str = f"{Fraction(y_val).limit_denominator()}"
                else:
                    y_str = f"{y_val:.6f}"
                table.add_row([Fore.YELLOW + x_str + Style.RESET_ALL, Fore.GREEN + y_str + Style.RESET_ALL])
            print(table)

        return {"x": self.x_points, "y": solution}


# Test cases and example usage
def test_solver():
    """
    Runs a suite of test cases to verify solver functionality.
    """
    # Test 1: Simple f(x) = 0 with Dirichlet boundary conditions
    def f1(x):
        return 0

    left_bc = BoundaryCondition(value=0)
    right_bc = BoundaryCondition(value=1)
    solver = FiniteDifferenceSolver(f1, a=0, b=1, n=10, left_bc=left_bc, right_bc=right_bc)
    result = solver.solve(verbose=True)
    expected = [i / 10 for i in range(11)]
    assert all(abs(r - e) < 1e-6 for r, e in zip(result["y"], expected)), "Test 1 failed."

    print(Fore.GREEN + "\nAll tests passed!")


if __name__ == "__main__":
    while True:
        print(Fore.CYAN + "\nFinite Difference Solver Menu:"
              "\n1. Run Test Cases"
              "\n2. Solve Custom Equation"
              "\n3. Exit")
        choice = input(Fore.YELLOW + "\nEnter your choice (1/2/3): ")

        if choice == '1':
            test_solver()
        elif choice == '2':
            try:
                a = float(input(Fore.YELLOW + "Enter the start of the interval (a): "))
                b = float(input(Fore.YELLOW + "Enter the end of the interval (b): "))
                n = int(input(Fore.YELLOW + "Enter the number of internal points (n): "))

                def custom_function(x):
                    return float(input(Fore.YELLOW + f"Enter the value of the function at x = {x}: "))

                left_bc_type = input(Fore.YELLOW + "Enter left boundary condition type (value/derivative): ").strip().lower()
                if left_bc_type not in ['value', 'derivative']:
                    raise ValueError("Invalid boundary condition type. Please enter either 'value' or 'derivative'.")
                if left_bc_type == 'value':
                    left_bc = BoundaryCondition(value=float(input(Fore.YELLOW + "Enter left boundary value: ")))
                elif left_bc_type == 'derivative':
                    left_bc = BoundaryCondition(derivative=float(input(Fore.YELLOW + "Enter left boundary derivative: ")))

                right_bc_type = input(Fore.YELLOW + "Enter right boundary condition type (value/derivative): ").strip().lower()
                if right_bc_type not in ['value', 'derivative']:
                    raise ValueError("Invalid boundary condition type. Please enter either 'value' or 'derivative'.")
                if right_bc_type == 'value':
                    right_bc = BoundaryCondition(value=float(input(Fore.YELLOW + "Enter right boundary value: ")))
                elif right_bc_type == 'derivative':
                    right_bc = BoundaryCondition(derivative=float(input(Fore.YELLOW + "Enter right boundary derivative: ")))

                # Show equation in a simplified format for confirmation
                equation_str = f"y''(x) = f(x)"
                boundary_str = f"Boundary Conditions: y(a) = {left_bc.value if left_bc.value is not None else 'N/A'}, y'(a) = {left_bc.derivative if left_bc.derivative is not None else 'N/A'}, " \
                                f"y(b) = {right_bc.value if right_bc.value is not None else 'N/A'}, y'(b) = {right_bc.derivative if right_bc.derivative is not None else 'N/A'}"
                print(Fore.MAGENTA + f"\nEquation: {equation_str}")
                print(Fore.MAGENTA + f"{boundary_str}")

                confirm = input(Fore.YELLOW + "\nIs this correct? (yes/no): ").strip().lower()
                if confirm != 'yes':
                    print(Fore.RED + "Operation canceled by the user.")
                    continue

                solver = FiniteDifferenceSolver(custom_function, a, b, n, left_bc, right_bc)
                solver.solve(verbose=True)

            except ValueError as e:
                print(Fore.RED + f"Error: {e}")
            except Exception as e:
                print(Fore.RED + f"Unexpected error: {e}")
        elif choice == '3':
            print(Fore.GREEN + "Exiting... Goodbye!")
            break
        else:
            print(Fore.RED + "Invalid choice. Please enter 1, 2, or 3.")
