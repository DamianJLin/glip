import warnings
import re
from lib.cyclic_list import CyclicList


class GaussDiagram():
    """
    Gauss diagram object, as used to compute the Tait graphs of a knot diagram from its Gauss code.
    """

    def __init__(self, gauss_code):
        gauss_words = re.findall('.*?[+-]', gauss_code)

        self.vertices = CyclicList()
        self.arcs = CyclicList()
        self.crossings = []

        # Construct GaussDiagram's vertices and crossings.
        crossing_number_to_crossing = {}
        for w in gauss_words:
            layer_input = w[0]
            crossing_number_input = w[1:-1]
            writhe_input = w[-1]

            if layer_input not in 'OU':
                raise ValueError("Error in gauss code. Words must start with 'O' or 'U'.")
            if not crossing_number_input.isdecimal():
                raise ValueError("Error is gauss code. Crossing indices must be decimal numbers.")
            if writhe_input not in '+-':
                raise ValueError("Error in gauss code. Words must end with '+' or '-'.")

            # TODO: Check that all labels appear twice and are in order.

            # Parse Gauss code, and adjust crossing numbers to start numbering from zero.
            layer = str(layer_input)
            crossing_number = int(crossing_number_input)
            if crossing_number == 0:
                warnings.warn(
                    "0 detected as crossing label in Gauss code. "
                    "Gauss codes should start labelling crossings from 1."
                )
            crossing_number -= 1
            writhe = {'+': 1, '-': -1}[writhe_input]

            # Create new Vertex and Crossing objects if necessary.
            v = self.GaussDiagramVertex(layer)
            if crossing_number in crossing_number_to_crossing:
                c = crossing_number_to_crossing[crossing_number]
            else:
                c = self.GaussDiagramCrossing(crossing_number, writhe)
                self.crossings.append(c)
                crossing_number_to_crossing[crossing_number] = c

            v.crossing = c
            if layer == 'O':
                c.vertex_over = v
            elif layer == 'U':
                c.vertex_under = v

            self.vertices.append(v)

        # Construct GaussDiagram's arcs.
        for i in range(len(self.vertices)):
            a = self.GaussDiagramArc()
            self.arcs.append(a)

            self.vertices[i].arc_next = a
            a.crossing_prev = self.vertices[i]
            self.vertices[i + 1].arc_prev = a
            a.crossing_next = self.vertices[i + 1]

    class GaussDiagramVertex():

        def __init__(self, layer):
            self.layer = layer

            self.crossing = None
            self.arc_next = None
            self.arc_prev = None

        def subsequent_arc(self, alignment):
            if alignment == 1:
                return self.arc_next
            elif alignment == -1:
                return self.arc_prev

        def jump_crossing(self):
            if self.layer == 'O':
                return self.crossing.vertex_under
            elif self.layer == 'U':
                return self.crossing.vertex_over

        def __repr__(self):
            layer_repr = self.layer
            if self.crossing is not None:
                crossing_number_repr = str(self.crossing.crossing_number)
                writhe_repr = {1: '+', -1: '-'}[self.crossing.writhe]
            else:
                crossing_number_repr = '_'
                writhe_repr = '_'
            return layer_repr + crossing_number_repr + writhe_repr

    class GaussDiagramCrossing():

        def __init__(self, crossing_number, writhe):
            self.crossing_number = crossing_number
            self.writhe = writhe

            self.vertex_over = None
            self.vertex_under = None

            self.adjacent_faces_green = [None, None]
            self.adjacent_faces_red = [None, None]

        def add_adjacent_face_green(self, face):
            if face is None:
                raise ValueError('Argument face cannot be None.')
            if self.adjacent_faces_green[0] is None:
                self.adjacent_faces_green[0] = face
            elif self.adjacent_faces_green[0] is not None and self.adjacent_faces_green[1] is None:
                self.adjacent_faces_green[1] = face
            elif self.adjacent_faces_green[0] is None and self.adjacent_faces_green[1] is None:
                raise AttributeError(f'{self}.adjacent_faces_green already has 2 elements.')

        def add_adjacent_face_red(self, face):
            if face is None:
                raise ValueError('Argument face cannot be None.')
            if self.adjacent_faces_red[0] is None:
                self.adjacent_faces_red[0] = face
            elif self.adjacent_faces_red[0] is not None and self.adjacent_faces_red[1] is None:
                self.adjacent_faces_red[1] = face
            elif self.adjacent_faces_red[0] is None and self.adjacent_faces_red[1] is None:
                raise AttributeError(f'{self}.adjacent_faces_red already has 2 elements.')

        def __repr__(self):
            crossing_number_repr = str(self.crossing_number)
            writhe_repr = {1: '+', -1: '-'}[self.writhe]
            return crossing_number_repr + writhe_repr

    class GaussDiagramArc():

        def __init__(self):
            self.crossing_next = None
            self.crossing_prev = None

        def subsequent_vertex(self, alignment):
            if alignment == 1:
                return self.crossing_next
            elif alignment == -1:
                return self.crossing_prev

        def __repr__(self):
            return f'({self.crossing_prev})--({self.crossing_next})'
