import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
from lib.cyclic_list import CyclicList
from lib.gauss_diagram import GaussDiagram


def gordon_litherland_green(gauss_code, verbose=False, very_verbose=False):
    """Compute the green Gordon Litherland matrix."""

    if very_verbose:
        verbose = True
    if verbose:
        sp.init_printing()

    gauss_diagram = GaussDiagram(gauss_code)

    # Check gauss code is alternating.
    if any(
        gauss_diagram.vertices[i].layer == gauss_diagram.vertices[i + 1].layer for i in
        range(len(gauss_diagram.vertices))
    ):
        raise ValueError(
            'Gauss code must be alternating, but {gauss_code} is not alternating. '
            'Find an alternating diagram (which every checkerboard colourable knot must have.)'
        )

    # === Step 1 ===

    # Traverse Gauss diagram for green Tait graph, computing edge list and crossing order, while
    # setting tokens.
    arcs_visited = set()
    face_counter = itertools.count()
    green_face_number_to_crossing_order = {}

    for i, arc in enumerate(gauss_diagram.arcs):
        # Choose initial turn of left for the green Tait graph.
        turn = (+1) * (-1) ** (i % 2)
        if arc not in arcs_visited:
            face_number = next(face_counter)
            crossing_order = CyclicList()
            arc_visit_order = CyclicList()
            alignment = 1
            vertex = gauss_diagram.vertices[i]
            vertex_0 = vertex

            # Traverse the face.
            not_finished = True
            while not_finished:
                crossing_order.append(vertex.crossing.crossing_number)
                # Move through an arc.
                arc = vertex.subsequent_arc(alignment)
                arcs_visited.add(arc)
                arc_visit_order.append(arc)
                vertex = arc.subsequent_vertex(alignment)
                # Move through a crossing.
                vertex.crossing.add_adjacent_face_green(face_number)
                vertex = vertex.jump_crossing()
                # Update alignment based on rule a = a * t * w * j where
                # t: 1/-1 for left/right
                # w: 1/-1 for pos/neg crossing sign
                # j: 1/-1 for over->under/under->over
                alignment *= turn * vertex.crossing.writhe * {'U': 1, 'O': -1}[vertex.layer]

                not_finished = vertex is not vertex_0 or alignment == -1

            # Reverse if necessary and set tokens.
            if turn == 1:
                for i, a in enumerate(arc_visit_order):
                    a.token = (face_number, i)
            elif turn == -1:
                crossing_order.reverse()
                arc_visit_order.reverse(style='reflect')
                for i, a in enumerate(arc_visit_order):
                    a.token = (face_number, i)

            green_face_number_to_crossing_order[face_number] = crossing_order

    green_tait_edge_list = {
        c.crossing_number: c.adjacent_faces_green for c in gauss_diagram.crossings
    }

    if verbose:
        print(f'green_tait_edge_list = \n{green_tait_edge_list}\n')
        print(f'green_face_number_to_crossing_order = \n{green_face_number_to_crossing_order}\n')

    # === Step 2 ===

    # Traverse Gauss diagram for red Tait graph, computing crossing order and boundary of red faces.
    arcs_visited = set()
    face_counter = itertools.count()
    red_face_number_to_crossing_order = {}
    red_face_number_to_token_order = {}

    for i, arc in enumerate(gauss_diagram.arcs):
        # Choose initial turn of left for the green Tait graph.
        turn = (-1) * (-1) ** (i % 2)
        if arc not in arcs_visited:
            face_number = next(face_counter)
            crossing_order = CyclicList()
            token_order = CyclicList()
            arc_visit_order = CyclicList()
            alignment = 1
            vertex = gauss_diagram.vertices[i]
            vertex_0 = vertex

            # Traverse the face.
            not_finished = True
            while not_finished:
                crossing_order.append(vertex.crossing.crossing_number)
                # Move through an arc.
                arc = vertex.subsequent_arc(alignment)
                token_order.append(arc.token)
                arcs_visited.add(arc)
                arc_visit_order.append(arc)
                vertex = arc.subsequent_vertex(alignment)
                # Move through a crossing.
                vertex.crossing.add_adjacent_face_red(face_number)
                vertex = vertex.jump_crossing()
                # Update alignment based on rule a = a * t * w * j where
                # t: 1/-1 for left/right
                # w: 1/-1 for pos/neg crossing sign
                # j: 1/-1 for over->under/under->over
                alignment *= turn * vertex.crossing.writhe * {'U': 1, 'O': -1}[vertex.layer]

                not_finished = vertex is not vertex_0 or alignment == -1

            # Reverse if necessary.
            if turn == 1:
                token_order.cycle(-1)
            if turn == -1:
                crossing_order.reverse()
                token_order.reverse()

            red_face_number_to_crossing_order[face_number] = crossing_order
            red_face_number_to_token_order[face_number] = token_order

    red_tait_edge_list = {c.crossing_number: c.adjacent_faces_red for c in gauss_diagram.crossings}

    if verbose:
        print(f'red_face_number_to_crossing_order = \n{red_face_number_to_crossing_order}\n')
        print(f'red_tait_edge_list = \n{red_tait_edge_list}\n')
        print(f'red_face_number_to_token_order = \n{red_face_number_to_token_order}\n')

    # === Step 3 ===

    # Compute the red 2-boundary operator and hence orient green edges.
    n_red_edges = len(red_tait_edge_list)
    n_red_faces = len(red_face_number_to_crossing_order)
    red_faces = dict(sorted(red_face_number_to_crossing_order.items())).values()
    red_2_boundary = np.zeros(
        shape=(n_red_edges, n_red_faces),
        dtype=int
    )

    green_face_number_to_orientation_order = {
        i: [None] * len(v) for i, v in green_face_number_to_crossing_order.items()
    }

    edges_seen = set()
    for face_number, face in enumerate(red_faces):
        for i, edge_number in enumerate(face):

            if edge_number not in edges_seen:
                # Compute red_2_boundary element.
                red_2_boundary[edge_number, face_number] += 1
                edges_seen.add(edge_number)
                # Set green edge orientation.
                green_face, green_ordinal = red_face_number_to_token_order[face_number][i]
                green_face_number_to_orientation_order[green_face][green_ordinal] = -1
            else:
                red_2_boundary[edge_number, face_number] -= 1
    # Set remaining green edge orientations (which have already been determined).
    for key in green_face_number_to_orientation_order:
        green_face_number_to_orientation_order[key] =\
            [1 if value is None else value
             for value in green_face_number_to_orientation_order[key]]

    # Use the red_2_boundary to calculate the kernel of the filling map (the inclusion of the red
    # spanning surface into the Carter surface).
    red_2_boundary = sp.Matrix(red_2_boundary)
    rank_ker = red_2_boundary.rank()
    basis_ker = red_2_boundary[:, :rank_ker]

    if verbose:
        print(f'red_2_boundary = \n{sp.pretty(red_2_boundary)}\n')
        print(
            f'green_face_number_to_orientation_order = \n'
            f'{green_face_number_to_orientation_order}\n'
        )
        print(f'basis_ker = \n{sp.pretty(basis_ker)}\n')

    # === Step 4 ===

    # Compute a basis for the homology of the green Tait graph by finding a maximally independent
    # set of the adjacency matrix.

    green_vertices = green_face_number_to_crossing_order  # Vertex/face duality.
    green_edges = sorted(green_tait_edge_list.keys())
    green_adjacency_matrix = np.zeros(
        shape=(len(green_vertices), len(green_edges)),
        dtype=int
    )

    # Compute the adjeacency matrix.
    for edge in green_edges:
        for vertex, cyclic_order in green_face_number_to_crossing_order.items():
            for i, c in enumerate(cyclic_order):
                if c == edge:
                    green_adjacency_matrix[vertex, edge] += \
                        green_face_number_to_orientation_order[vertex][i]

    green_adjacency_matrix = sp.Matrix(green_adjacency_matrix)
    basis_green_tait_homology = sp.Matrix.hstack(
        *green_adjacency_matrix.nullspace()
    )

    if verbose:
        print(f'green_adjacency_matrix = \n{sp.pretty(green_adjacency_matrix)}\n')
        print(f'basis_green_tait_homology = \n{sp.pretty(basis_green_tait_homology)}\n')

    # === Step 5 ===

    # Extend basis_ker to a basis for the green Tait homology.
    rank = sp.shape(basis_green_tait_homology)[1]
    basis_extended = np.zeros(
        shape=(len(green_edges), rank),
        dtype=int
    )
    basis_extended = sp.Matrix(basis_extended)
    rank_asymmetric = rank - rank_ker
    basis_extended[:, rank_asymmetric:] = basis_ker  # Set the symmetric part on the end.

    i = 0  # The next column of basis_green_tait_homology to try.
    j = 0  # The next empty column of basis_extended
    while basis_extended.rank() < rank:
        old_rank = basis_extended.rank()
        basis_extended[:, j] = basis_green_tait_homology[:, i]
        new_rank = basis_extended.rank()
        if old_rank == new_rank:
            # Convert back.
            basis_extended[:, j] = sp.zeros(len(green_edges), 1)
            i += 1
        elif old_rank < new_rank:
            i += 1
            j += 1

    if verbose:
        print(f'basis_extended = \n{sp.pretty(basis_extended)}\n')

    # === Step 6 ===

    # Compute the symmetric term of the form.

    form_symmetric = basis_extended.transpose() @ basis_extended

    if verbose:
        print(f'form_symmetric = \n{sp.pretty(form_symmetric)}\n')

    # === Step 7 ===

    # Compute the antisymmetric term of the form.
    form_antisymmetric = sp.Matrix.zeros(rank)

    if very_verbose:
        print('Local linking numbers will be printed. For reference:')
        print(f'green_face_number_to_crossing_order = \n{green_face_number_to_crossing_order}\n')
        print(
            f'green_face_number_to_orientation_order = \n'
            f'{green_face_number_to_orientation_order}\n'
        )

    # Compute the antisymmetric form.
    for i, j in itertools.combinations(range(rank), 2):

        linking_number = 0
        cycle_a = basis_extended[:, i]
        cycle_b = basis_extended[:, j]
        edges_in_a = set()
        edges_in_b = set()
        for e in green_edges:
            if cycle_a[e] != 0:
                edges_in_a.add(e)
            if cycle_b[e] != 0:
                edges_in_b.add(e)

        # Compute reduced version of green_face_number_to_crossing_order and
        # green_face_number_to_orientation_order.
        green_crossing_dict_reduced = {}
        green_orientation_dict_reduced = {}
        for vertex, crossing_order in green_face_number_to_crossing_order.items():
            green_crossing_dict_reduced[vertex] = CyclicList()
            green_orientation_dict_reduced[vertex] = CyclicList()
            for k, e in enumerate(crossing_order):
                if e in edges_in_a or e in edges_in_b:
                    green_crossing_dict_reduced[vertex].append(e)
                    green_orientation_dict_reduced[vertex].append(
                        green_face_number_to_orientation_order[vertex][k]
                    )

        # Compute local linking number at each vertex.
        for vertex in green_crossing_dict_reduced:
            crossing_order = green_crossing_dict_reduced[vertex]
            orientation_order = green_orientation_dict_reduced[vertex]
            #  Switch to convention 1 for out, -1 for in.
            orientation_order = CyclicList(*[-o for o in orientation_order])

            # Split into cases of degree 2, 3 and 4.
            degree = len(green_orientation_dict_reduced[vertex])

            if very_verbose:
                print(f'\ta: {edges_in_a}\t{sp.pretty(cycle_a.transpose())}')
                print(f'\tb: {edges_in_b}\t{sp.pretty(cycle_b.transpose())}')
                print(f'\tvertex: {vertex}')
                print(f'\tdegree: {degree}')
                print('\tcrossing_order / orientation_order:')
                print('\t', crossing_order)
                print('\t', orientation_order)

            # TODO: Properly handle cases where degree is greater than 4 using recursion.

            if degree == 0:
                local_linking_number = 0

            elif degree == 2:
                local_linking_number = 0

            elif degree == 3:
                k = 0
                while not crossing_order[0] in set.intersection(edges_in_a, edges_in_b):
                    assert k <= 2, 'k Assertion Failed'
                    crossing_order.cycle()
                    orientation_order.cycle()
                    k += 1

                if very_verbose:
                    print(f'\t{k = }')

                local_linking_number = Fraction(1, 2)

                if crossing_order[1] in edges_in_a:
                    local_linking_number *= 1
                else:
                    local_linking_number *= -1

                # Account for orientation of edges.
                local_linking_number *= orientation_order[1] * orientation_order[2]
                # Find orientation of the cycles with respect to edges, that is
                # account for negative coefficients in cycle_a, cycle_b.
                cycle_1 = cycle_a if crossing_order[1] in edges_in_a else cycle_b
                cycle_2 = cycle_a if crossing_order[2] in edges_in_a else cycle_b
                local_linking_number *= cycle_1[crossing_order[1]] * cycle_2[crossing_order[2]]

            elif degree == 4:

                # Confirm both cycles have some intersection.
                if not set.intersection(edges_in_a, crossing_order)\
                        or not set.intersection(edges_in_b, crossing_order):
                    continue

                k = 0
                while not (
                    crossing_order[0] in edges_in_a
                    and
                    orientation_order[0] * cycle_a[crossing_order[0]] == 1
                ):
                    assert k <= 3, 'k Assertion Failed'
                    crossing_order.cycle()
                    orientation_order.cycle()
                    k += 1

                if very_verbose:
                    print(f'\t{k = }')

                # Account for orientation of edges
                if crossing_order[1] in edges_in_a or crossing_order[-1] in edges_in_a:
                    local_linking_number = 0
                elif crossing_order[1] in edges_in_b:
                    if orientation_order[1] == 1:
                        local_linking_number = 1
                    elif orientation_order[1] == -1:
                        local_linking_number = -1
                # Find orientation of the cycles with respect to edges, that is
                # account for negative coefficients in cycle_a, cycle_b.
                cycle_0 = cycle_a if crossing_order[0] in edges_in_a else cycle_b
                cycle_1 = cycle_a if crossing_order[1] in edges_in_a else cycle_b
                local_linking_number *= cycle_0[crossing_order[0]] * cycle_1[crossing_order[1]]

            else:
                raise NotImplementedError('degrees odd or 3 is not yet implemented')

            if very_verbose:
                print(f'\t{local_linking_number = }\n')

            linking_number += local_linking_number

        form_antisymmetric[i, j] = linking_number
        form_antisymmetric[j, i] = -linking_number

    if verbose:
        print(f'form_antisymmetric = \n{sp.pretty(form_antisymmetric)}\n')

    # === Step 8 ===

    form = form_symmetric + form_antisymmetric

    if verbose:
        print(f'form = \n{sp.pretty(form)}\n')

    return form
