import itertools
import functools
import numpy as np
import sympy as sp
from lib.cyclic_list import CyclicList
from lib.gauss_diagram import GaussDiagram
from lib.basis_ke import basis_ker_ext
import lib.error as er
import lib.linking as lk


def gordon_litherland(gauss_code, symmetric=False, verbose=False, very_verbose=False):
    """
    Compute the matrices of the Gordon Litherland pairing.
    If symmetric is True, compute instead the Symmetric Gordon Litherland Pairing.
    """

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
        raise er.GaussCodeNotAlternatingError(
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

    # Compute the red 2-boundary operator and hence orient green/red edges.
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
    red_face_number_to_orientation_order = {
        i: [None] * len(v) for i, v in red_face_number_to_crossing_order.items()
    }

    edges_seen = set()
    for face_number, face in enumerate(red_faces):
        for i, edge_number in enumerate(face):

            if edge_number not in edges_seen:
                # Compute red_2_boundary element.
                red_2_boundary[edge_number, face_number] += 1
                edges_seen.add(edge_number)
                # Set green/red edge orientation.
                red_face_number_to_orientation_order[face_number][i] = -1
                green_face, green_ordinal = red_face_number_to_token_order[face_number][i]
                green_face_number_to_orientation_order[green_face][green_ordinal] = -1
            else:
                red_2_boundary[edge_number, face_number] -= 1
    # Set remaining green/red edge orientations (which have already been determined).
    for key in green_face_number_to_orientation_order:
        green_face_number_to_orientation_order[key] =\
            [1 if value is None else value
             for value in green_face_number_to_orientation_order[key]]
    for key in red_face_number_to_orientation_order:
        red_face_number_to_orientation_order[key] =\
            [1 if value is None else value
             for value in red_face_number_to_orientation_order[key]]

    # Use green tait graph orientation to also compute green_2_boundary.
    n_green_edges = len(green_tait_edge_list)
    n_green_faces = len(green_face_number_to_crossing_order)
    green_faces = dict(sorted(green_face_number_to_crossing_order.items())).values()
    green_2_boundary = np.zeros(
        shape=(n_green_edges, n_green_faces),
        dtype=int
    )
    for face_number, face in enumerate(green_faces):
        for i, edge_number in enumerate(face):
            green_2_boundary[edge_number, face_number] +=\
                -green_face_number_to_orientation_order[face_number][i]

    # === Step 4 ===

    # Compute a basis for the homology of the green Tait graph by finding a maximally independent
    # set of the adjacency matrix.

    green_basis = basis_ker_ext(
        green_face_number_to_crossing_order,
        green_face_number_to_orientation_order,
        green_tait_edge_list,
        red_2_boundary
    )

    red_basis = basis_ker_ext(
        red_face_number_to_crossing_order,
        red_face_number_to_orientation_order,
        red_tait_edge_list,
        green_2_boundary
    )

    # === Step 6 ===

    # Compute the symmetric term of the form.

    def form_symmetric(basis):
        return basis.transpose() @ basis

    green_form_symmetric = form_symmetric(green_basis)
    red_form_symmetric = form_symmetric(red_basis)

    if verbose:
        print(f'green_form_symmetric = \n{sp.pretty(green_form_symmetric)}\n')
        print(f'red_form_symmetric = \n{sp.pretty(red_form_symmetric)}\n')

    if symmetric:
        return (green_form_symmetric, red_form_symmetric)

    # === Step 7 ===

    def form_antisymmetric(basis, crossing_orders, orientation_orders):

        rank = sp.shape(basis)[1]
        form_antisymmetric = sp.zeros(rank)

        linking_number = functools.partial(
            lk.linking_number,
            crossing_orders=crossing_orders,
            orientation_orders=orientation_orders,
            very_verbose=very_verbose
        )

        if very_verbose:
            print('Local linking numbers will be printed. For reference:')
            print(f'crossing_orders = \n{crossing_orders}\n')
            print(f'orientation_orders = \n{orientation_orders}\n')

        for i, j in itertools.combinations(range(rank), 2):

            cycle_a = basis[:, i]
            cycle_b = basis[:, j]
            linking_number_ij = linking_number(cycle_a=cycle_a, cycle_b=cycle_b)

            form_antisymmetric[i, j] = linking_number_ij
            form_antisymmetric[j, i] = -linking_number_ij

        return form_antisymmetric

    green_form_antisymmetric = form_antisymmetric(
        green_basis,
        green_face_number_to_crossing_order,
        green_face_number_to_orientation_order
    )
    red_form_antisymmetric = form_antisymmetric(
        red_basis,
        red_face_number_to_crossing_order,
        red_face_number_to_orientation_order
    )

    if verbose:
        print(f'green_form_antisymmetric = \n{sp.pretty(green_form_antisymmetric)}\n')
        print(f'red_form_antisymmetric = \n{sp.pretty(red_form_antisymmetric)}\n')

    # === Step 8 ===

    green_form = green_form_symmetric + green_form_antisymmetric
    red_form = red_form_symmetric + red_form_antisymmetric

    if verbose:
        print(f'red_form = \n{sp.pretty(red_form)}\n')
        print(f'green_form = \n{sp.pretty(green_form)}\n')

    return (red_form, green_form)
