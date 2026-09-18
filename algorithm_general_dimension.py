"""
Algorithm to compute the lattice diameter and a lattice diameter segment of a rational polytope.

Function:
  Input: A is a Sage matrix over QQ, and b is a Sage vector over QQ.
  * lattice_diameter_ILP(A, b):
    returns the lattice diameter and one lattice diameter segment, returns (-Infinity, None) if P has no lattice points, and raise ValueError if P contains lattice points but is unbounded.
  
  Alternatively, for a list of vertices run
  * lattice_diameter_ILP_by_vertices(vertices)

Algorithm:
For fixed integer s, decide whether there are integer vectors x, y, u with 
    y - x = s*u,     x, y in P,     and u_i >= 1
for at least one coordinate i. Feasibility means that the lattice diameter is
at least s. Binary search then finds the
largest feasible s.

The implementation requires SageMath. The integer linear programs are solved using SageMath’s PPL backend, which uses exact rational arithmetic.
"""

from sage.all import Infinity, Polyhedron, QQ, matrix, vector
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException



# ==================================================
# HELPER FUNCTIONS
# ==================================================

def _add_membership_constraints(model, var, A, b):
    """Add A*x <= b to the model."""
    for i in range(A.nrows()):
        model.add_constraint(sum(A[i, j] * var[j] for j in range(A.ncols())) <= b[i])


def _any_lattice_point(A, b):
    """Return one integer point z satisfying A*z <= b, or None if none exists."""
    model = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    z = model.new_variable(integer=True, nonnegative=False)
    _add_membership_constraints(model, z, A, b)
    model.set_objective(0)

    try:
        model.solve()
    except MIPSolverException:
        return None

    values = model.get_values(z)
    return tuple(int(values[j]) for j in range(A.ncols()))

def _feasible_for_coordinate(A, b, s, coordinate):
    """
    Find integer vectors x, y, u with A*x <= b, A*y <= b, y - x = s*u,
    and u[coordinate] >= 1. Return (x, y, u), or None if none exist.
    """
    dimension = A.ncols()
    model = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    x = model.new_variable(integer=True, nonnegative=False)
    y = model.new_variable(integer=True, nonnegative=False)
    u = model.new_variable(integer=True, nonnegative=False)

    _add_membership_constraints(model, x, A, b)
    _add_membership_constraints(model, y, A, b)

    for j in range(dimension):
        model.add_constraint(y[j] - x[j] == s * u[j])
    model.add_constraint(u[coordinate] >= 1)
    model.set_objective(0)

    try:
        model.solve()
    except MIPSolverException:
        return None

    x_values = model.get_values(x)
    y_values = model.get_values(y)
    u_values = model.get_values(u)
    return (
        tuple(int(x_values[j]) for j in range(dimension)),
        tuple(int(y_values[j]) for j in range(dimension)),
        tuple(int(u_values[j]) for j in range(dimension)),
    )

def _feasible(A, b, s):
    """
    s will be a positive integer.
    Return an integer witness (x, y, u) that the lattice diameter is at least s,
    or None if it is smaller than s. 
    """
    for coordinate in range(A.ncols()):
        solution = _feasible_for_coordinate(A, b, s, coordinate)
        if solution is not None:
            return solution
    return None

def _check_unbounded(A):
    """
    Checks whether recession cone {r : Ar <= 0} contains nonzero vector, i.e, if
        - rank(A) < dimension, or
        - Ar <= 0 and r nonzero, equivalently at least one entry of Ar is negative, 
            equivalently (after scaling r) sum_i (Ar)_i <= -1
    """
    dimension = A.ncols()
    if A.rank() < dimension:
        return True

    # find nonzero r with A*r <= 0, via Ar <=0 and sum_i (Ar)_i <= -1
    model = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    r = model.new_variable(real=True, nonnegative=False)
    rows = [
        sum(A[i, j] * r[j] for j in range(dimension))
        for i in range(A.nrows())
    ]
    for row in rows:
        model.add_constraint(row <= 0)
    model.add_constraint(sum(rows) <= -1)
    model.set_objective(0) #feasibility problem
    try:
        model.solve()
    except MIPSolverException as exc:
        message = str(exc).lower()
        if "infeasible" in message or "no feasible solution" in message:
            return False
        raise
    return True


# ==================================================
# MAIN
# - lattice_diameter_ILP
# - lattice_diameter_ILP_by_vertices
# ==================================================

def lattice_diameter_ILP(A, b):
    """
    Input: A is a Sage matrix over QQ, and b is a Sage vector over QQ.
    Returns the lattice diameter and a lattice diameter segment of P = {x : Ax <= b}.
    If there are no lattice points, return (-Infinity, None). If it's unbounded, raise ValueError.
    """
    if len(b) != A.nrows() or A.ncols() == 0:
        raise ValueError("A and or b are the wrong dimension.")

    lattice_point = _any_lattice_point(A, b) #check if P contains a lattice point
    if lattice_point is None:
        return -Infinity, None
    
    if _check_unbounded(A):
        raise ValueError("Expected a bounded polytope")
    
    #check whether lattice diameter >= 1
    solution = _feasible(A, b, 1)
    if solution is None:
        return 0, (lattice_point, lattice_point)

    lower = 1
    best_solution = solution
    test_value = 2

    #finds largest k with lattice diameter >= 2^k, set upper = 2^{k+1}-1, lower = 2^k
    while True:
        solution = _feasible(A, b, test_value)
        if solution is None:
            upper = test_value - 1
            break
        lower = test_value
        best_solution = solution
        test_value *= 2

    #do binary search till upper = lower 
    while lower < upper:
        middle = (lower + upper + 1) // 2
        solution = _feasible(A, b, middle)
        if solution is None:
            upper = middle - 1
        else:
            lower = middle
            best_solution = solution

    x, y, _ = best_solution
    return lower, (x, y)

def lattice_diameter_ILP_by_vertices(vertices):
    """
    Input: vertices of a rational polytope.
    Returns the lattice diameter and a lattice-diameter segment.

    """
    P = Polyhedron(vertices=vertices, base_ring=QQ)

    inequalities = P.inequalities_list()
    equations = P.equations_list()

    rows = [[-coefficient for coefficient in inequality[1:]] for inequality in inequalities]
    bounds = [inequality[0] for inequality in inequalities]

    for equation in equations:
        constant = equation[0]
        coefficients = list(equation[1:])
        rows.append(coefficients)
        bounds.append(-constant)
        rows.append([-coefficient for coefficient in coefficients])
        bounds.append(constant)

    A = matrix(QQ, rows)
    b = vector(QQ, bounds)

    return lattice_diameter_ILP(A, b)