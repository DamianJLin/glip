import re
import itertools


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

        self.jump = None
        self.writhe = writhe

    def __repr__(self):
        height_str = {1: 'O', -1: 'U'}[self.height]
        index_str = str(self.index)
        writhe_str = {1: '+', -1: '-'}[self.writhe]
        return height_str + index_str + writhe_str

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
                v.jump = u
                u.jump = v
                connect_to[v.index - 1] = None
        if any(i is not None for i in connect_to):
            raise ValueError(
                "gauss_code contains a crossing index with no pair."
            )

        self.n = max(v.index for v in self.vertices)
        assert self.n * 2 == len(self.vertices)

    def tait_graphs(self):
        """
        Takes ChordDiagram out returns tuple of two tait graphs for the
        corresponding knot.
        """
        graphs = []

        for initial_turn in (-1, 1):
            strand_visited = [None] * self.n * 2
            edge_dict = {i + 1: [] for i in range(self.n)}
            counter = itertools.count()

            def strand_index_from_vertex(vertex: ChordVertex, alignment):
                if alignment == 1:
                    return v.pos_in_gc
                elif alignment == -1:
                    return (v.pos_in_gc - 1) % (self.n * 2)

            for i, s in enumerate(strand_visited):
                turn = initial_turn * (-1) ** (i % 2)
                if s is None:
                    face = next(counter)
                    # Set starting alignment of 1.
                    a = 1
                    v = self.vertices[i]
                    while strand_visited[strand_index_from_vertex(v, a)]\
                            is None:
                        # Set strand to be visited by current face.
                        strand_visited[strand_index_from_vertex(v, a)] = face
                        # Tranverse strand.
                        v = v.traverse(a)
                        # Set crossing as visited by current face.
                        edge_dict[v.index].append(face)
                        # Jump crossing.
                        v = v.jump
                        # Update alignment.
                        a = a * turn * v.height * v.writhe
            assert all(len(faces) == 2 for faces in edge_dict.values())
            graphs.append(tuple(edge_dict.values()))

        return tuple(graphs)


def tait_graph_edge_lists(gauss_code):
    return ChordDiagram(gauss_code).tait_graphs()
