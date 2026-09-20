# lattice-diameters

- `algorithm_2d_edge_opposite_vertex.py`: Code to compute the lattice diameter of a lattice polygon P and a lattice diameter segment for every lattice diameter direction of P.

- `algorithm_general_dimension.py`: Code to compute the lattice diameter of any rational polytope given either by inequalities with rational entries or by its vertices.

- `plotting.py`: Code for plotting polytopes and their lattice diameter segments; used in `examples.ipynb`.

- `examples.ipynb`: Demonstrates the algorithms on a lattice polygon and a three-dimensional lattice polytope.

## Requirements

This project requires [SageMath](https://doc.sagemath.org/html/en/installation/).

## Installation

Clone the repository and enter its directory:

```bash
git clone https://github.com/aebrose/lattice-diameters.git
cd lattice-diameters
```

Install SageMath according to its [installation guide](https://doc.sagemath.org/html/en/installation/).

## Example usage

Create a Python file in the `lattice-diameters` directory. For example, save the following as `my_example.py`:

```python
from algorithm_2d_edge_opposite_vertex import lattice_diameter_2d, lattice_diameter_2d_all_directions

vertices = [(0, 1), (1, -2), (-1, -2), (-3, 0), (-2, 1)]

ld, ld_seg = lattice_diameter_2d(vertices)
_, ld_segs_dir = lattice_diameter_2d_all_directions(vertices)

print("Lattice diameter:", ld)
print("A lattice diameter segment:", ld_seg)
print("Lattice diameter segments, at least one for every lattice diameter direction:", ld_segs_dir)
```

or, for a higher-dimensional example:

```python
from sage.all import QQ, matrix, vector
from algorithm_general_dimension import lattice_diameter_ILP, lattice_diameter_ILP_by_vertices

A = matrix(QQ, [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
b = vector(QQ, [2, 0, 2, 0, 2, 0])

ld, ld_seg = lattice_diameter_ILP(A, b)

print("Lattice diameter:", ld)
print("A lattice diameter segment:", ld_seg)

#Or, if the rational polytope is given by vertices:

vertices = [(0,0,0), (2,0,1), (0,QQ(3)/2,0), (0,0,QQ(5)/3)]

ld, ld_seg = lattice_diameter_ILP_by_vertices(vertices)

print("Lattice diameter:", ld)
print("A lattice diameter segment:", ld_seg)
```

Then run the file with SageMath:

```bash
sage -python my_example.py
```

See `examples.ipynb` for two- and three-dimensional examples with plots.
