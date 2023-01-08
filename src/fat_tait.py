import re
import itertools


class FatEdgeList():
    """
    Edge list that keeps track of cyclic ordering of edges.
    """

    def __init__(self, edge_list, cord_dict):
        self.edge_list = edge_list
        self.cord_dict = cord_dict


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

    def fat_tait_graphs(self):
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
            # Create empty dictionary to keep track of cyclic order.
            cyclic_ord_dict = {}
            # Counter to assign labels to faces.
            counter = itertools.count()

            # Function to get the 'next' strand based on the current vertex in
            # the Gauss diagram traversal.
            def next_strand_idx(vertex: ChordVertex, alignment):
                if alignment == 1:
                    return v.pos_in_gc
                elif alignment == -1:
                    return (v.pos_in_gc - 1) % (self.n * 2)

            # This loops over strands and if the strand is not yet visited,
            # finds the corresponding face.
            for i, visited in enumerate(strand_visited):
                turn = initial_turn * (-1) ** (i % 2)
                if not visited:
                    face = next(counter)
                    # Set starting alignment of 1.
                    a = 1
                    v = self.vertices[i]
                    cyclic_ord = []
                    while not strand_visited[next_strand_idx(v, a)]:
                        # Set strand as visited.
                        strand_visited[next_strand_idx(v, a)] = True
                        # Append the current edge (vertex of chord diagram)
                        # cyclic order.
                        cyclic_ord.append(v.index)
                        # Caluclate jump variable (+1/-1 for over to
                        # under/under to over).
                        jump = v.height
                        # Move around the chord diagram.
                        v = v.traverse(a)
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
                        cyclic_ord.reverse()
                    cyclic_ord_dict[face] = cyclic_ord
            assert all(len(faces) == 2 for faces in edge_list.values())
            graphs.append(
                FatEdgeList(tuple(edge_list.values()), cyclic_ord_dict)
            )

        return tuple(graphs)


def get_fat_tait_graphs(gauss_code):
    return ChordDiagram(gauss_code).fat_tait_graphs()
