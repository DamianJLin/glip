import fat_tait as ft
import numpy as np
import scipy as scp
import sympy as sp

np.set_printoptions(sign=' ')

# Algorithm as follows
# Step 1: Compute basis for ker i* by taking a maximal linearly independent set
# of faces of the opposite colour in the knot diagram.
# Step 2: Extend this basis to a basis for the 1st homology of the Tait graph.
# Step 3: Calculate the symmetric part of the GL-form on the extended basis.
# This is equivalent to the inner product (Gram matrix) of the lattice of
# integer flows.


def gl_pairing(gauss_code, swap_graphs=False):

    # === Setup ===

    g, h = ft.get_fat_tait_graphs(gauss_code)

    if swap_graphs:
        g, h = h, g

    faces = list(h.cord_dict.values())
    edges = g.edge_list

    # Edge labels inside faces start from 1 but we need them to start from 0 to
    # correspond to the correct edge in edges.
    for i, face in enumerate(faces):
        faces[i] = [edge - 1 for edge in face]
    print(f'g: {g.edge_list}, {g.cord_dict},\nh: {h.edge_list}, {h.cord_dict}')

    # Calculate the number of n-cells in the Carter surface (n_edges,
    # n_vertices also agree with the Tait graph).
    vertices = set()
    for edge in edges:
        vertices |= set(edge)
    n_vertices = len(vertices)

    n_edges = len(edges)

    n_faces = len(faces)

    # === Step 1 ===

    # Construct the 2-boundary operator on the faces of the Carter surface with
    # columns faces written in the edge basis.
    carter_boundary_2 = np.zeros(
        shape=(n_edges, n_faces),
        dtype=int
    )
    edge_seen = [False] * n_edges
    for j, face in enumerate(faces):

        for i, edge in enumerate(face):
            if not edge_seen[edge]:
                carter_boundary_2[edge, j] += 1
                edge_seen[edge] = True
            else:
                carter_boundary_2[edge, j] -= 1

    print(f'carter_boundary_2 = \n{carter_boundary_2}')

    rank_ker = np.linalg.matrix_rank(carter_boundary_2)
    ker_subbasis = carter_boundary_2[:, :rank_ker]
    print(f'ker_subbasis = \n{ker_subbasis}')

    # === Step 2 ===

    # Construct the 1-boundary operator of the Tait graph with columns edges
    # written in the vertex basis.
    tait_boundary_1 = np.zeros(
        shape=(n_vertices, n_edges),
        dtype=int
    )
    for j, edge in enumerate(edges):
        # final - initial
        tait_boundary_1[edge[1], j] += 1
        tait_boundary_1[edge[0], j] -= 1
    print(f'tait_boundary_1 = \n{tait_boundary_1}')

    tait_homology_basis = sp.Matrix(tait_boundary_1).nullspace()
    # A transposition is necessary converting between numpy as sympy.
    tait_homology_basis = np.asarray(
        tait_homology_basis, dtype=int).squeeze().transpose()

    print(f'tait_homology_basis = \n{tait_homology_basis}')

    # Extend ker_subbasis to a basis for the homology of the Tait graph.

    rank = tait_homology_basis.shape[1]
    basis = np.zeros(
        shape=(n_edges, rank),
        dtype=int
    )
    rank_asym = rank - rank_ker
    basis[:, rank_asym:] = ker_subbasis

    i = 0  # The next column of tait_homology_basis to fill.
    j = 0  # The next empty column of basis.
    while np.linalg.matrix_rank(basis) < rank:
        # Compute candidate for basis extension.
        ext_cand = basis.copy()
        ext_cand[:, j] = tait_homology_basis[:, i]
        if np.linalg.matrix_rank(basis) < np.linalg.matrix_rank(ext_cand):
            basis = ext_cand.copy()
            j += 1
        i += 1

    print(
        f'Basis matrix for the antisymmetric subspace '
        f'(of dim {rank_asym}): \n{basis[:, :rank_asym]}'
    )
    print(f'Extension to H_1(T) (of dim {rank}): \n{basis[:, :rank]}')

    # === Step 3 ===

    gl_sym = basis.transpose() @ basis
    print(f'gl_sym = \n{gl_sym}')

    # === Step 4 ===

    gl_asym = np.zeros_like(gl_sym)

    # Implement the traversal algorithm to augment the fat Tait graph to
    # include orientation information.

    fat_lists = list(g.cord_dict.values())

    # Edge labels inside faces start from 1 but we need them to start from 0 to
    # correspond to the correct edge in edges.
    for i, fat_list in enumerate(fat_lists):
        fat_lists[i] = [edge - 1 for edge in fat_list]
    print(fat_lists)

    def find_other(pos):
        vertex, index = pos
        assert index < len(fat_lists[vertex])
        target = fat_lists[vertex][index]
        for j, fat_list in enumerate(fat_lists):
            for i, edge in enumerate(fat_list):
                if (j, i) != (vertex, index) and fat_lists[j][i] == target:
                    return j, i
        ValueError('Failed to find other.')

    def clockwise(pos):
        vertex, index = pos
        return vertex, (index - 1) % len(fat_lists[vertex])


# Testing. TODO: Remove.
if __name__ == '__main__':
    # 5_2429: O1-U2-O3+U1-O2-U4-O5-U3+O4-U5-
    # 5_2428: O1-U2-O3-U1-O2-U4+O5+U3-O4+U5+
    gl_pairing('O1-U2-O3+U1-O2-U4-O5-U3+O4-U5-')
