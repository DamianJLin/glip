import sympy as sp
import src.fat_tait as ft


def convert_to_finite_field(matrix, order):
    m, n = matrix.shape
    field = sp.GF(order)
    rows = []
    for i in range(m):
        row = []
        for j in range(n):
            row.append(field(matrix[i, j]))
        rows.append(row)
    return sp.polys.matrices.DomainMatrix(rows, matrix.shape, field)


def virtual_genus(gauss_code):

    g, h = ft.get_fat_tait_graphs(gauss_code)

    # Compute dim H_1(X, Z_2)
    vertices = set()

    for edge in g.edge_list:
        vertices |= set(edge)
    m = len(vertices)
    n = len(g.edge_list)

    boundary_1 = sp.zeros(m, n)
    for j, edge in enumerate(g.edge_list):
        # final - initial
        boundary_1[edge[0], j] += 1
        boundary_1[edge[1], j] -= 1

    boundary_1 = convert_to_finite_field(boundary_1, 2)
    dim_homology_G = n - boundary_1.rank()

    # Compute dim ker_i
    # We count the dimension of the span of the disks of opposite colour in the
    # Carter surface (which we know generate the kernel), in order to find the
    # dimension of that kernel.
    faces = list(h.cord_dict.values())
    edges = set()
    for face in faces:
        edges |= set(face)
    m = len(edges)
    n = len(faces)

    boundary_2 = sp.zeros(m, n)
    for j, face in enumerate(faces):
        for i in face:
            # We are in Z_2 do we don't need to worry about orientation.
            # i - 1 since face labels start from 1.
            boundary_2[i - 1, j] += 1

    boundary_2 = convert_to_finite_field(boundary_2, 2)
    dim_ker_i = boundary_2.rank()

    # We should now have the dimension of the Carter surface.
    dim_homology_carter = dim_homology_G - dim_ker_i
    carter_genus = dim_homology_carter // 2
    return carter_genus
