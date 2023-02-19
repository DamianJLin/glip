import warnings
import re
from lib.cyclic_list import CyclicList


class GaussDiagramVertex():

    def __init__(self, layer):
        self.layer = layer

        self.crossing = None
        self.arc_next = None
        self.arc_prev = None

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

    def __repr__(self):
        crossing_number_repr = str(self.crossing_number)
        writhe_repr = {1: '+', -1: '-'}[self.writhe]
        return crossing_number_repr + writhe_repr


class GaussDiagramArc():

    def __init__(self):
        self.crossing_next = None
        self.crossing_prev = None

    def __repr__(self):
        return str(self.crossing_prev) + '--' + str(self.crossing_next)


class GaussDiagram():

    def __init__(self, gauss_code):
        gauss_words = re.findall('.*?[+-]', gauss_code)

        self.vertices = CyclicList()

        # Construct GaussDiagram's vertices and crossings.
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

            v = GaussDiagramVertex(layer)
            c = GaussDiagramCrossing(crossing_number, writhe)

            v.crossing = c
            if layer == 'O':
                c.vertex_over = v
            elif layer == 'U':
                c.vertex_under = v

            self.vertices.append(v)

        # Construct GaussDiagram's arcs.
        for i in range(len(self.vertices)):
            a = GaussDiagramArc()

            self.vertices[i].arc_next = a
            a.crossing_prev = self.vertices[i]
            self.vertices[i + 1].arc_prev = a
            a.crossing_next = self.vertices[i + 1]
