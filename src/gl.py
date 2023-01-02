import fat_tait as ft
import collections
import numpy as np

# Algorithm as follows
# Step 1: Compute basis for ker i by taking a maximal linearly independent set
# of faces of the opposite colour in the knot diagram.
# Step 2: Extend this basis.


def gl_pairing(gauss_code):

    g, h = ft.get_fat_tait_graphs(gauss_code)
    g, h = h, g

    faces = list(h.cord_dict.values())
    edges = g.edge_list

    # Edge labels inside faces start from 1 but we need them to start from 0 to
    # correspond to the correct edge in edges.
    for i, face in enumerate(faces):
        faces[i] = [edge - 1 for edge in face]
    print(f'g: {g.edge_list}, {g.cord_dict},\nh: {h.edge_list}, {h.cord_dict}')

    # Calculate dimensions of vector spaces
    vertices = set()
    for edge in edges:
        vertices |= set(edge)
    n_vertices = len(vertices)

    n_edges = len(edges)

    n_faces = len(faces)

    # Construct the 1-boundary operator with columns edges written in the
    # vertex basis.
    boundary_1 = np.zeros(
        shape=(n_vertices, n_edges),
        dtype=int
    )
    for j, edge in enumerate(edges):
        # final - initial
        boundary_1[edge[1], j] += 1
        boundary_1[edge[0], j] -= 1
    print(f'boundary_1 = \n{boundary_1}')

    # Construct the 2-boundary operator with columns faces written in the
    # edge basis.
    boundary_2 = np.zeros(
        shape=(n_edges, n_faces),
        dtype=int
    )
    for j, face in enumerate(faces):
        # To make sure we choose correct orientation for edges, made sure
        # d^2 = 0; i.e. boundary_1 of boundary_2 is 0.
        # TODO: This is dodgy. Find a way to orient self-loops in g's Tait
        # graph to avoid this hack.
        curr = np.zeros(n_edges)

        # Remove items with any duplicates but retain order.
        face = [k for k, v in collections.Counter(face).items() if v == 1]

        for i, edge in enumerate(face):
            nxt = np.zeros(n_edges)
            nxt[edge] = 1
            d_next = boundary_1 @ nxt
            d_curr = boundary_1 @ curr
            prod = - np.dot(d_curr, d_next)
            if prod == 0:
                inc = 1
            elif np.abs(prod) == 1:
                inc = prod
            elif np.abs(prod) == 2:
                inc = np.sign(prod)

            curr += inc * nxt

        boundary_2[:, j] = curr
    print(f'boundary_2 = \n{boundary_2}')

    print(f'd^2 = \n{boundary_1 @ boundary_2}')


# Testing. TODO: Remove.
if __name__ == '__main__':
    gl_pairing('O1-U2-O3+U1-O2-U4-O5-U3+O4-U5-')
