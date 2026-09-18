"""
Algorithm to compute lattice diameter segments of lattice polygons and all lattice diameter directions by iterating through all edge-opposite vertex triangles. 

Functions:
    Input: P_vertices is a finite list of integral vectors whose convex hull defines the lattice polygon.
  * lattice_diameter_2d(P_vertices): 
        returns the lattice diameter and one lattice diameter segment.
  * lattice_diameter_2d_all_directions(P_vertices):
        returns the lattice diameter, and at least one lattice diameter segment for every lattice diameter direction

Algorithm:
1. Construct the polygon, sort its vertices and edges in cyclic order.
2. For every edge e with outward normal a, find its opposite vertex (or opposite vertices) if there is an edge parallel to e. That is, v is an opposite vertex to e, if for an outward normal vector a of e, v minimizes the inner product <a,-> over all points in P.
3. For all edge-opposite vertex triangles T = conv(e,v):
    - Find up to three local lattice diameters of T at v by minimizing <a,-> over lattice points in T (but not v). This is done by solving (up to three) integer linear programs.
4. A lattice diameter segment is a local lattice diameter segment of some T at v, with the maximum lattice length. All lattice diameter directions are realized by a local lattice diameter segment of some T at v.

The implementation requires SageMath. The integer linear programs are solved using SageMath’s PPL backend, which uses exact rational arithmetic.
"""

from functools import cmp_to_key
from sage.all import Polyhedron, QQ, ZZ, gcd
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

# ==================================================
# PREPROCESSING: Helper functions, and defines a class for `LatticePolygon` and `Edge` storing necessary data for algorithm.
# ==================================================

def dot(a, v):
    return a[0] * v[0] + a[1] * v[1]

def cyclic_order(vertices):
    """ Sorts vertices in counter-clockwise cyclic order based on direction from centroid."""
    n = len(vertices)
    center = tuple(QQ(sum(v[i] for v in vertices)) / n for i in range(2))

    def compare_angle(p, q, center):
        """
        Input: two points p,q and a center point (in the interior of the polygon)
        Output:
            * -1 if p comes before q in counter clockwise ordering, starting at (1,0). For example, if p in top halfspace, q in bottom halfspace, or if in same halfspace and cross product (p-c) x (q-c) > 0
            * +1 if q comes before p
        """
        dir_p = (p[0] - center[0], p[1] - center[1])
        dir_q = (q[0] - center[0], q[1] - center[1])
        p_half = 0 if (dir_p[1] > 0 or (dir_p[1] == 0 and dir_p[0] > 0)) else 1
        q_half = 0 if (dir_q[1] > 0 or (dir_q[1] == 0 and dir_q[0] > 0)) else 1
        if p_half != q_half:
            return -1 if p_half < q_half else 1
        cross = dir_p[0] * dir_q[1] - dir_p[1] * dir_q[0]
        if cross > 0:
            return -1
        else: 
            return 1
    
    compare = lambda p, q: compare_angle(p, q, center)
    return tuple(sorted(vertices, key=cmp_to_key(compare)))


class LatticePolygon:
    """
    Define an instance of LatticePolygon by LatticePolygon(P), for a list of vertices.
    The class preprocess the lattice polygon with useful data for the algorithm.
        - stores vertices and edges in counter-clockwise order
        - precomputes opposite vertex/vertices for each edge
        - stores edges as an instance of the class Edge, with attributes .vertices, .a (outward normal) and .offset (<v,a> for v in e).
    """

    class Edge:
        """ Lets us access outward normal vector a and offset by e.a and e.offset."""
        def __init__(self, vertices, a, offset):
            self.vertices = tuple(vertices)
            self.a = tuple(a)
            self.offset = int(offset)

    def __init__(self, vertices): 
        self.polyhedron = Polyhedron(vertices=[tuple(v) for v in vertices], base_ring=ZZ)
        if self.polyhedron.dim() != 2:
            raise ValueError(f"expected a 2-dimensional polygon, got dim {self.polyhedron.dim()}")
        V = tuple(tuple(int(c) for c in v) for v in self.polyhedron.vertices_list())
        self.vertices_ccw = cyclic_order(V)
        self._edges = self._make_edges()
        self._opposites = self._precompute_opposites()

    def _make_edges(self):
        """
        Makes each edge `e` an instance of `Edge' which stores its vertices, outward normal vector and offset.
        """
        edges = []
        n = len(self.vertices_ccw)
        for i, p in enumerate(self.vertices_ccw):
            q = self.vertices_ccw[(i + 1) % n]
            a = (q[1] - p[1], p[0] - q[0]) #outward normal vector, since q comes after p in ccw order
            edges.append(self.Edge((p, q), a, dot(a,p)))
        return tuple(edges)

    def _precompute_opposites(self):
        """Compute all opposite vertices to edges. By walking through vertices in cyclic ordering it never revisits a vertex and has O(n) operations."""
        n = len(self.vertices_ccw)
        first_a = self._edges[0].a
        first_values = [dot(first_a,v) for v in self.vertices_ccw]
        first_minimum = min(first_values)

        # Finds opposite vertex to edge e0, if there are two it picks the second one (in ccw order).
        j = next(k for k in range(n)
                 if first_values[k] == first_minimum and first_values[(k + 1) % n] > first_minimum)
        
        stop = j + n
        opposites = {}

        for e in self._edges:
            #while loop finds index j (in cyclic order) that minimizes <e.a,->, if there are two then it picks second one
            while j < stop:
                current_val = dot(e.a, self.vertices_ccw[j % n])
                next_val = dot(e.a, self.vertices_ccw[(j + 1) % n])
                if next_val > current_val: #current is minimizer, i.e opposite vertex to the edge
                    break
                j += 1

            #while loop breaks when we find minimum value
            min_val = dot(e.a, self.vertices_ccw[j % n])
            previous_val = dot(e.a, self.vertices_ccw[(j - 1) % n])
            if previous_val == min_val: #two opposite vertices
                minimizing_vertices = (self.vertices_ccw[(j - 1) % n], self.vertices_ccw[j % n])
            else: #one opposite vertex 
                minimizing_vertices = (self.vertices_ccw[j % n],)
            opposites[e] = (min_val, minimizing_vertices)

        return opposites

    def edges(self):
        return self._edges

    def opposite_vertices(self, edge):
        return self._opposites[edge]

# ==================================================
# ILP & Helpers to call ILP
# ==================================================

def rows_from_vertices(vertices):
    """
    Returns inequality rows for the H-description defining conv(vertices) used for ILP solver. It makes facet normals primitive, this is not necessary but experiments show that on average it speeds up the ILP solver."""
    polyhedron = Polyhedron(vertices=[tuple(v) for v in vertices], base_ring=ZZ)
    rows = []
    for inequality in polyhedron.inequalities_list():
        b, a1, a2 = (int(value) for value in inequality)
        divisor = int(gcd(a1, a2)) #makes ILP a bit faster 
        rows.append((a1 // divisor, a2 // divisor, b // divisor))
    return tuple(rows)

def halfspace_row(a, rhs):
    """The row for the ILP solver encoding  <a, x> >= rhs,  i.e.  a1*x + a2*y - rhs >= 0."""
    return (int(a[0]), int(a[1]), -int(rhs))


def solve_ilp(rows, objective):
    """
    Minimize <objective, x> over the lattice points x in Z^2 satisfying `rows`.
    Input:  rows: a list of tuples (a1, a2, b), meaning  a1*x + a2*y + b >= 0.
    objective : (c1, c2) integer coefficients.
    Output: (optimal_value, optimal_solution)
    """
    p = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    x = p.new_variable(integer=True, nonnegative=False)
    for (a1, a2, b) in rows:
        p.add_constraint(a1 * x[0] + a2 * x[1] + b >= 0)
    p.set_objective(objective[0] * x[0] + objective[1] * x[1])
    try:
        opt_val = p.solve()
    except MIPSolverException as error:
        raise RuntimeError("The triangle ILP unexpectedly failed.") from error 
    opt_sol = p.get_values(x)
    return opt_val, (int(opt_sol[0]), int(opt_sol[1]))


# ==================================================
# MAIN
# - lattice_diameter_2d:
#   computes lattice diameter, and one lattice diameter segment
# - lattice_diameter_2d_all_directions:
#   computes lattice diameter, and at least one lattice diameter segment per lattice diameter direction
# ==================================================

def lattice_diameter_2d(P_vertices):
    """
    Input: list of integral vectors defining the polygon. 
    Returns: (ld, ld_seg)
        where ld is the lattice diameter and ld_seg is a lattice diameter segment of P.
    """
    P = LatticePolygon(P_vertices) #makes instance that precomputes edge-opposite vertices and useful notation
    ld = 0
    for e in P.edges():
        opp_value, opp_verts = P.opposite_vertices(e)
        for v in opp_verts:
            rows = list(rows_from_vertices(list(e.vertices) + [v])) + [halfspace_row(e.a, opp_value + 1)] #H-rep of triangle and extra ILP inequality
            aw, w = solve_ilp(rows, objective=e.a)
            lattice_length = int((e.offset - opp_value) // (aw - opp_value))
            if lattice_length > ld:
                ld = lattice_length
                ld_seg = (v, (v[0] + lattice_length * (w[0] - v[0]),
                                v[1] + lattice_length * (w[1] - v[1])))
    return ld, ld_seg


def other_points_at_same_height(w, a, rows):
    """
    Helper for lattice_diameter_2d_all_directions:
    Returns up to two other lattice points z, with <a,z> = <a,w>, satisfying constraints in `rows`.
    """
    g = int(gcd(a[0], a[1]))
    step = (a[1] // g, -a[0] // g) #primitive edge direction 
    points = []
    #at most 3 local lattice diameter segments, hence at most two neighbors to w
    for step_size in [-1,1,-2,2]: 
        if len(points) < 2:
            z = (w[0] + step_size * step[0], w[1] + step_size * step[1])
            if all(a1 * z[0] + a2 * z[1] + b >= 0 for a1, a2, b in rows):
                points.append(z)
    return points

def lattice_diameter_2d_all_directions(P_vertices):
    """
    Input: list of integral vectors defining the polygon P
    Returns: (ld, ld_segs)
        where ld is the lattice diameter of P, and ld_segs contains at least one lattice diameter segment per lattice diameter direction (but possibly more)
    """
    P = LatticePolygon(P_vertices)
    ld = 0
    ld_segs = []

    for e in P.edges():
        opp_value, opp_verts = P.opposite_vertices(e)

        for v in opp_verts:
            num_ld_at_v = 0
            maybe_more = False
            triangle_rows = list(rows_from_vertices(list(e.vertices) + [v]))
            rows = triangle_rows + [halfspace_row(e.a, opp_value + 1)]
            aw, w = solve_ilp(rows, objective=e.a)
            lattice_length = int((e.offset - opp_value) // (aw - opp_value))

            if lattice_length > ld:
                maybe_more = True
                num_ld_at_v = 1
                ld = lattice_length
                ld_segs = [(v, (v[0] + lattice_length * (w[0] - v[0]),
                                v[1] + lattice_length * (w[1] - v[1])))]
            elif lattice_length == ld:
                num_ld_at_v = 1
                maybe_more = True
                ld_segs.append((v, (v[0] + lattice_length * (w[0] - v[0]),
                                    v[1] + lattice_length * (w[1] - v[1]))))
                
            while maybe_more and num_ld_at_v < 3:
                points = other_points_at_same_height(w, e.a, triangle_rows)
                ld_segs.extend((v, (v[0] + lattice_length * (p[0] - v[0]),
                                    v[1] + lattice_length * (p[1] - v[1]))) for p in points)
                num_ld_at_v += len(points)
                if aw + 1 <= e.offset and num_ld_at_v < 3:
                    rows = list(triangle_rows) + [halfspace_row(e.a, aw + 1)]
                    aw, w = solve_ilp(rows, objective=e.a)
                    ll = int((e.offset - opp_value) // (aw - opp_value))
                    if ll == lattice_length:
                        num_ld_at_v += 1
                        ld_segs.append((v, (v[0] + lattice_length * (w[0] - v[0]),
                                            v[1] + lattice_length * (w[1] - v[1]))))
                    else: break
                else: break
    return ld, ld_segs
