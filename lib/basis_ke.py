import sympy as sp
import numpy as np


def basis_ker_ext(crossing_orders, orientation_orders, edge_list, opp_2_boundary):
    """
    Compute a basis for the 1st homology of the Tait graph with colour cl
    with the first k columns being a basis for the asymmetric subspace
    and the last n - k columns being a basis for the kernel of the filling map.
    """

    # Compute a basis for the kernel of the filling map.
    opp_2_boundary = sp.Matrix(opp_2_boundary)
    rank_ker = opp_2_boundary.rank()
    basis_ker = opp_2_boundary[:, :rank_ker]

    # Compute adjacency matrix.
    vertices = crossing_orders  # By vertex/face duality.
    edges = sorted(edge_list.keys())
    adjacency_matrix = np.zeros(
        shape=(len(vertices), len(edges)),
        dtype=int
    )

    for edge in edges:
        for vertex, crossing_order in crossing_orders.items():
            for i, c in enumerate(crossing_order):
                if c == edge:
                    adjacency_matrix[vertex, edge] += orientation_orders[vertex][i]

    adjacency_matrix = sp.Matrix(adjacency_matrix)
    basis_homology = sp.Matrix.hstack(
        *adjacency_matrix.nullspace()
    )

    # Extend basis_ker to basis for homology.
    rank = sp.shape(basis_homology)[1]
    basis_extended = np.zeros(
        shape=(len(edges), rank),
        dtype=int
    )
    basis_extended = sp.Matrix(basis_extended)
    rank_asymmetric = rank - rank_ker
    basis_extended[:, rank_asymmetric:] = basis_ker

    i = 0  # The next column of basis_homology to try.
    j = 0  # The next empty column of basis_extended.
    while basis_extended.rank() < rank:
        old_rank = basis_extended.rank()
        basis_extended[:, j] = basis_homology[:, i]
        new_rank = basis_extended.rank()
        if old_rank == new_rank:
            # Convert back
            basis_extended[:, j] = sp.zeros(len(edges), 1)
            i += 1
        elif old_rank < new_rank:
            i += 1
            j += 1

    return basis_extended
