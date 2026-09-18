"""Plot a two- or three-dimensional polytope with one or many of its lattice diameter segments."""

from sage.all import Polyhedron, ceil, floor, gcd, line, line3d, point, point3d


def lattice_points_on_segment(segment):
    p, q = segment
    differences = [int(q[i] - p[i]) for i in range(len(p))]
    number_of_steps = 0
    for difference in differences:
        number_of_steps = int(gcd(number_of_steps, abs(difference)))
    if number_of_steps == 0:
        return [tuple(p)]
    step = [difference // number_of_steps for difference in differences]
    return [tuple(p[i] + k * step[i] for i in range(len(p)))
            for k in range(number_of_steps + 1)]


def plot_lattice_diameter(ld_segs, vertices=None, A=None, b=None):
    """
    Input: a list of lattice diameter segments, and either a list of vertices or a pair A, b
        defining the polytope by A*x <= b.
    Returns: a Sage graphics object showing the polytope and lattice diameter segment.
    """
    if vertices is not None:
        P = Polyhedron(vertices=vertices)
    elif A is not None and b is not None:
        inequalities = [[b[i]] + [-value for value in A.row(i)] for i in range(A.nrows())]
        P = Polyhedron(ieqs=inequalities)

    if P.dim() not in (2, 3) or P.dim() != P.ambient_dim():
        raise ValueError("Expected a 2- or 3-dimensional polytope")
    if not P.is_compact():
        raise ValueError("expected a bounded polytope")

    segment_colors = ["orange", "red", "green", "purple"]

    if P.dim() == 2:
        polytope_vertices = P.vertices_list()
        x_min = int(floor(min(v[0] for v in polytope_vertices)))
        x_max = int(ceil(max(v[0] for v in polytope_vertices)))
        y_min = int(floor(min(v[1] for v in polytope_vertices)))
        y_max = int(ceil(max(v[1] for v in polytope_vertices)))
        lattice_points = [(x, y) for x in range(x_min, x_max+1)
                           for y in range(y_min, y_max+1)]

        picture = P.plot(
            point=False, #vertices
            line={"color": "black", "thickness": 1}, #edges
            polygon={"color": "lightgray", "alpha": 0.35},
            axes = False
        )
        picture += point(lattice_points, color="lightgray", size=12)
        for i, s in enumerate(ld_segs):
            color = segment_colors[i % len(segment_colors)]
            picture += line(s, color=color, thickness=3)
            picture += point(lattice_points_on_segment(s), color=color, size=70)
        return picture

    polytope_vertices = P.vertices_list()
    x_min = int(floor(min(v[0] for v in polytope_vertices)))
    x_max = int(ceil(max(v[0] for v in polytope_vertices)))
    y_min = int(floor(min(v[1] for v in polytope_vertices)))
    y_max = int(ceil(max(v[1] for v in polytope_vertices)))
    z_min = int(floor(min(v[2] for v in polytope_vertices)))
    z_max = int(ceil(max(v[2] for v in polytope_vertices)))
    lattice_points = [(x, y, z) for x in range(x_min, x_max + 1)
                       for y in range(y_min, y_max + 1)
                       for z in range(z_min, z_max + 1)]

    picture = P.plot(
        point=False,
        line={"color": "black", "thickness": 1},
        polygon={"color": "lightblue", "opacity": 0.35},
        frame=False
    )
    picture += point3d(lattice_points, color="gray", size=15, opacity=0.9)
    for i, s in enumerate(ld_segs):
        color = segment_colors[i % len(segment_colors)]
        picture += line3d(s, color=color, thickness=2)
        picture += point3d(lattice_points_on_segment(s), color=color, size=15)
    return picture
