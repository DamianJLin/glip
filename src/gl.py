import numpy as np
import sympy as sp
import itertools
import re

np.set_printoptions(sign=' ')


# This is really messy right now - as happens when you continually
# expand an algorithm as you go.

# Algorithm as follows:

# Step 0: Modified version of the Fat tait graph algorithm that keeps track
# of information later used to orient non-trivial homology generators in
# step 4.

# Step 1: Compute basis for ker i* by taking a maximal linearly independent set
# of faces of the opposite colour in the knot diagram.

# Step 2: Extend this basis to a basis for the 1st homology of the Tait graph.

# Step 3: Calculate the symmetric part of the GL-form on the extended basis.
# This is equivalent to the inner product (Gram matrix) of the lattice of
# integer flows.

# Step 4: Compute the asymmetric part of the GL-form.


# === Step 0 ===

class FatEdgeListPlus():
    """
    Edge list that keeps track of cyclic ordering of edges.
    """

    def __init__(self, edge_list, crossing_cycl_ord, alignment_ord,
                 overstrand_pair_pos):
        self.edge_list = edge_list
        self.crossing_cycl_ord = crossing_cycl_ord
        self.alignment_ord = alignment_ord
        self.overstrand_pair_pos = overstrand_pair_pos


class ChordVertex():
    """
    Vertex object on a ChordDiagram
    """

    def __init__(self, index, height, writhe):
        self.index = index
        self.height = height
        self.pos_in_gc = None

        self.traverse_next = None
        self.traverse_last = None

        self.cross = None
        self.writhe = writhe

    def __repr__(self):
        height_str = {1: 'O', -1: 'U'}[self.height]
        index_str = str(self.index)
        writhe_str = {1: '+', -1: '-'}[self.writhe]
        return height_str + index_str + writhe_str

    # Move one edge around the chord diagram.
    def traverse(self, alignment):
        if alignment == 1:
            return self.traverse_next
        elif alignment == -1:
            return self.traverse_last
        else:
            raise ValueError()


class ChordDiagram():
    """
    Object representing a chord diagram (cyclic graph on words in gauss code
    with edges between vertices representing the same crossing) imbued with
    some additional methods and structure.
    """

    def __init__(self, gauss_code: str):
        gauss_code = re.findall('.*?[+-]', gauss_code)
        self.vertices = []

        # Create ChordVertex objects.
        for i, s in enumerate(gauss_code):
            index = int(s[1:-1])
            height = {'O': 1, 'U': -1}[s[0]]
            writhe = {'+': 1, '-': -1}[s[-1]]
            v = ChordVertex(index, height, writhe)
            v.pos_in_gc = i
            self.vertices.append(v)

        # Create traversal edges and set pos_in_gc.
        for u, v in zip(self.vertices, self.vertices[1:] + self.vertices[:1]):
            u.traverse_next = v
            v.traverse_last = u

        # Create crossing edges.
        connect_to = [None] * len(self.vertices)
        for i, v in enumerate(self.vertices):
            if connect_to[v.index - 1] is None:
                connect_to[v.index - 1] = i
            else:
                u = self.vertices[connect_to[v.index - 1]]
                v.cross = u
                u.cross = v
                connect_to[v.index - 1] = None
        if any(i is not None for i in connect_to):
            raise ValueError(
                "gauss_code contains a crossing index with no pair."
            )

        self.n = max(v.index for v in self.vertices)
        assert self.n * 2 == len(self.vertices)

    def fat_tait_graphs_plus(self):
        """
        Takes ChordDiagram and returns tuple of two tait graphs for the
        corresponding knot.
        """

        # Knots are checkerboard colourable in Carter surface iff they are
        # alternating. Otherwise Tait graph is not defined. We check.
        init_height = self.vertices[0].height
        for i, v in enumerate(self.vertices):
            if not init_height == (-1) ** (i) * v.height:
                raise ValueError('Gauss code is not alternating.')

        # For each Tait graph, we create an array of all the strands that
        # have been visited. For each strand not yet visited, we explore
        # the face that the strand belongs to, keeping track of the
        # connections between faces diagonally across crossings, and the
        # cyclic order of the crossings around each face. As we visit, we
        # mark all discovered strands as visited. When all strands are
        # visited we are done. This is done via the chord diagram.

        # There will be two fat Tait graphs. Keep list to add them to.
        graphs = []

        # The turn direction dictates which checkerboard colour we are looking
        # at, and hence which Tait graph.
        for initial_turn in (1, -1):

            # Create array of which strands have been visited.
            strand_visited = [False] * self.n * 2
            # Create empty dictionary associating each crossing to the faces it
            # connects. Will be filled as we traverse.
            edge_list = {i + 1: [] for i in range(self.n)}
            # Create empty dictionary to keep track of cyclic order. Augmented
            # to also maintain alignment.
            cyclic_crossings_dict = {}
            cyclic_alignment_dict = {}
            # Counter to assign labels to faces.
            counter = itertools.count()

            # Overstrand association map (the aumgmentation to the algorithm).
            # Will map crossing number to a choice of two options for
            # (face, cyclic_id) based on alignment.
            stub_pos_at_overstrand = []
            for i in range(self.n):
                stub_pos_at_overstrand.append({})

            # Function to get the 'next' strand based on the current vertex in
            # the Gauss diagram traversal.
            def next_strand_idx(vertex: ChordVertex, alignment):
                if alignment == 1:
                    return v.pos_in_gc
                elif alignment == -1:
                    return (v.pos_in_gc - 1) % (self.n * 2)

            # This loops over strands and if the strand is not yet visited,
            # finds the corresponding face.
            ctr = 0
            for i, visited in enumerate(strand_visited):
                ctr += 1
                turn = initial_turn * (-1) ** (i % 2)
                if not visited:
                    face = next(counter)
                    # Set starting alignment of 1.
                    a = 1
                    v = self.vertices[i]
                    cyclic_crossings = []
                    cyclic_alignment = []
                    while not strand_visited[next_strand_idx(v, a)]:
                        # Set next strand as visited.
                        strand_visited[next_strand_idx(v, a)] = True
                        # Set stub_pos_at_overstrand.

                        if v.height == 1:  # Overstrand
                            # TODO: See why this is bugged and does not
                            # complete for the second Tait graph (though it
                            # won't affect the algorithm unless swap_graphs is
                            # set to True).
                            # print(f'{v.index = }, {a = }')
                            stub_pos_at_overstrand[v.index - 1][a] = \
                                (face, len(cyclic_crossings))
                            # print(f'{stub_pos_at_overstrand = }')
                        # Append the current edge (vertex of chord diagram)
                        # cyclic order.
                        if initial_turn == -1 and ctr == 2:
                            print('/', v.index - 1)
                            print(f'here: {turn * a}')
                        cyclic_crossings.append(v.index)
                        # TODO: Check that * -v.writhe works in all cases
                        cyclic_alignment.append(a * turn * -v.writhe)
                        # Move around the chord diagram.
                        v = v.traverse(a)
                        # Caluclate jump variable (+1/-1 for over to
                        # under/under to over).
                        jump = v.height
                        # Set crossing as visited by current face.
                        edge_list[v.index].append(face)
                        # Cross the chord diagram.
                        v = v.cross
                        # Update alignment.
                        a = a * turn * jump * v.writhe

                    # If turn direction was left then the cyclic order was
                    # correct. Otherwise it's the reversed and we need
                    # to correct it.
                    if turn == -1:
                        cyclic_crossings.reverse()
                        cyclic_alignment.reverse()
                        cyclic_alignment = cyclic_alignment[1:] + \
                            cyclic_alignment[:1]
                        print([c - 1 for c in cyclic_crossings],
                              cyclic_alignment)
                    cyclic_crossings_dict[face] = cyclic_crossings
                    cyclic_alignment_dict[face] = cyclic_alignment

            assert all(len(faces) == 2 for faces in edge_list.values())

            graphs.append(
                FatEdgeListPlus(
                    tuple(edge_list.values()),
                    cyclic_crossings_dict,
                    cyclic_alignment_dict,
                    stub_pos_at_overstrand
                )
            )

        return tuple(graphs)


def get_fat_tait_graphs_plus(gauss_code):
    return ChordDiagram(gauss_code).fat_tait_graphs_plus()


def gl_pairing(gauss_code, swap_graphs=False):

    # === Setup ===

    g, h = get_fat_tait_graphs_plus(gauss_code)

    print(g.overstrand_pair_pos)

    if swap_graphs:
        g, h = h, g

    faces = list(h.crossing_cycl_ord.values())
    edges = g.edge_list

    # Edge labels inside faces start from 1 but we need them to start from 0 to
    # correspond to the correct edge in edges.
    for i, face in enumerate(faces):
        faces[i] = [edge - 1 for edge in face]
    print(
        f'g: {g.edge_list}, {g.crossing_cycl_ord},\n'
        f'h: {h.edge_list}, {h.crossing_cycl_ord}'
    )

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
        extension_candidate = basis.copy()
        extension_candidate[:, j] = tait_homology_basis[:, i]
        if np.linalg.matrix_rank(basis) \
                < np.linalg.matrix_rank(extension_candidate):
            basis = extension_candidate.copy()
            j += 1
        i += 1

    print(
        f'Basis matrix for the antisymmetric subspace '
        f'(of dim {rank_asym}): \n{basis[:, :rank_asym]}'
    )
    print(f'Extension to H_1(T) (of dim {rank}): \n{basis[:, :rank]}')

    # === Step 3 ===

    form_symmetric = basis.transpose() @ basis
    print(f'form_symmetric = \n{form_symmetric}')

    # === Step 4 ===

    form_asymmetric = np.zeros_like(form_symmetric)

    # Implement the traversal algorithm to augment the fat Tait graph to
    # include orientation information.

    # This is what we want to be able to read off orientations from.
    fat_lists = list(g.crossing_cycl_ord.values())
    # Edge labels inside faces start from 1 but we need them to start from 0 to
    # correspond to the correct edge in edges.
    for i, fat_list in enumerate(fat_lists):
        fat_lists[i] = [edge - 1 for edge in fat_list]
    print(f'{fat_lists = }')

    # This is what we need to compute that.
    alignment_dict = h.alignment_ord
    position_association = g.overstrand_pair_pos

    # Compute in/out orientations around each vertex.
    fat_orienations = []
    for ls in fat_lists:
        fat_orienations.append([None] * len(ls))
    edge_seen = [False] * n_edges
    print(faces)
    print(f'{alignment_dict = }')
    print(f'{position_association = }')
    for j, face in enumerate(faces):
        for i, edge in enumerate(face):
            print(f'{edge}')
            if not edge_seen[edge]:
                al = alignment_dict[j][i]
                v, cyclic_id = position_association[edge][al]
                fat_orienations[v][cyclic_id] = -1
                print(f'{al = }')
                print(f'{fat_orienations=}')
                edge_seen[edge] = True

    # Now that the "in's" have been recorded, the rest are "out's".
    for fo in fat_orienations:
        for i, e in enumerate(fo):
            if e is None:
                fo[i] = 1

    print(f'{fat_orienations = }')

    # Now that we have the orientations we can actally calculate the asymmetric
    # form.
    # Loop over all unique pairs of nontrivial homology generating cycles.
    # TODO: Convince yourself these only have 1 "in and 1 "out" per vertex.
    for j in range(rank_asym):
        for i in range(j):
            cycle_a = basis[:, i]
            cycle_b = basis[:, j]
            print(cycle_a, cycle_b)
            # TODO: Stop hardcoding this.
            for k, fat_list in enumerate(fat_lists):
                e1 = None
                e2 = None
                m_0 = None
                for m, c in enumerate(fat_list):
                    if cycle_a[c] != 0:
                        print(f'{m = }')
                        e1_cycle = 'a'
                        e1 = - cycle_a[c] * fat_orienations[k][m]
                        m_0 = m
                        break
                    elif cycle_b[c] != 0:
                        print(f'{m = }')
                        e1_cycle = 'b'
                        e1 = - cycle_b[c] * fat_orienations[k][m]
                        m_0 = m
                        break
                for m, c in enumerate(fat_list):
                    if m == m_0:
                        continue
                    if cycle_a[c] != 0:
                        print(f'{m = }')
                        e2_cycle = 'a'
                        e2 = - cycle_a[c] * fat_orienations[k][m]
                        break
                    elif cycle_b[c] != 0:
                        print(f'{m = }')
                        e2_cycle = 'b'
                        e2 = - cycle_b[c] * fat_orienations[k][m]
                        break
                if e1_cycle == e2_cycle:
                    break
                coefficient = e1 * e2 * {'a': 1, 'b': -1}[e1_cycle]
                print(coefficient)
                form_asymmetric[i, j] = 1 * coefficient
                form_asymmetric[j, i] = -1 * coefficient

    print(f'form_asymmetric =\n{form_asymmetric}')

    # === Terminate ===

    form = form_symmetric + form_asymmetric
    print(
        f'The Gordon-Litherland form of knot with Gauss code {gauss_code} is:'
    )
    print(form)

    def inv(mat):
        mat = sp.Matrix(mat)
        return (mat.transpose() @ mat.inverse_CH()).trace()
    print(
        f'Alleged invariant, tr(M^TM^-1), for testing: \n{inv(form)}'
    )

    print('Invariant of Hans\' gl form for 5_2429:')
    print(inv(np.asarray(
        [
            [1, 1, 0, -1, 0],
            [-1, 1, -1, 0, 0],
            [0, -1, 2, 0, -1],
            [-1, 0, 0, 2, -1],
            [0, 0, -1, -1, 3]
        ]
    )))


# Testing. TODO: Remove.
if __name__ == '__main__':
    # 5_2429: O1-U2-O3+U1-O2-U4-O5-U3+O4-U5-
    # 5_2428: O1-U2-O3-U1-O2-U4+O5+U3-O4+U5+
    gl_pairing('O1-U2-O3+U1-O2-U4-O5-U3+O4-U5-')
