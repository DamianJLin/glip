from lib.cyclic_list import CyclicList
from fractions import Fraction
import sympy as sp


class Bracket:
    def __init__(self, cycle, direction):
        if cycle not in ('a', 'b'):
            raise ValueError("Argument cycle must have value 'a' or 'b'.")
        if direction not in (1, 0, -1):
            raise ValueError("Argument cycle must have value -1 or 1, or 0 for the empty bracket.")
        self.cycle = cycle
        self.direction = direction

    def __repr__(self):
        rep = self.cycle + {1: '>', -1: '<', 0: '-'}[self.direction]
        return rep


def linking_number(crossing_orders, orientation_orders, cycle_a, cycle_b, very_verbose):
    edges = range(sp.shape(cycle_a)[0])
    edgeset_a = set()
    edgeset_b = set()
    for e in edges:
        if cycle_a[e] != 0:
            edgeset_a.add(e)
        if cycle_b[e] != 0:
            edgeset_b.add(e)

    # Filter out edges not in a or b.
    crossing_orders_reduced = {}
    orientation_orders_reduced = {}

    for vertex, crossing_order in crossing_orders.items():
        crossing_orders_reduced[vertex] = CyclicList()
        orientation_orders_reduced[vertex] = CyclicList()

        for k, e in enumerate(crossing_order):
            edgeset_ab = set.union(edgeset_a, edgeset_b)
            if e in edgeset_ab:
                crossing_orders_reduced[vertex].append(e)
                orientation_orders_reduced[vertex].append(
                    orientation_orders[vertex][k]
                )

    def local_linking_number(vertex):
        crossing_order = crossing_orders_reduced[vertex]
        orientation_order = orientation_orders_reduced[vertex]
        # Switch orientation to -1/1 for in/out.
        orientation_order = CyclicList(*[-o for o in orientation_order])

        bracketing_a_priority = []
        bracketing_b_priority = []

        for e, o in zip(crossing_order, orientation_order):
            bracket_a = Bracket('a', cycle_a[e] * o)
            bracket_b = Bracket('b', cycle_b[e] * o)

            if bracket_a.direction != 0:
                bracketing_a_priority.append(bracket_a)
            if bracket_b.direction != 0:
                bracketing_a_priority.append(bracket_b)

            if bracket_b.direction != 0:
                bracketing_b_priority.append(bracket_b)
            if bracket_a.direction != 0:
                bracketing_b_priority.append(bracket_a)

        # Compute linking numberings given by each bracketing and average.
        level_a = 0
        local_linking_number_a = 0
        local_linking_number_b = 0
        for bk in bracketing_a_priority:
            if bk.cycle == 'a':
                level_a += bk.direction
            if bk.cycle == 'b':
                local_linking_number_a += bk.direction * level_a
        level_a = 0
        for bk in bracketing_b_priority:
            if bk.cycle == 'a':
                level_a += bk.direction
            if bk.cycle == 'b':
                local_linking_number_b += bk.direction * level_a

        local_linking_number = Fraction(local_linking_number_a + local_linking_number_b, 2)

        if very_verbose:
            print(f'\ta: {edgeset_a}\t{sp.pretty(cycle_a.transpose())}')
            print(f'\tb: {edgeset_b}\t{sp.pretty(cycle_b.transpose())}')
            print(f'\tvertex: {vertex}')
            print(f'\tbracketing_a_priority: {bracketing_a_priority}')
            print('\tcrossing_order / orientation_order:')
            print('\t', crossing_order)
            print('\t', orientation_order)
            print('\t', f'local_linking_number: {local_linking_number}')
            print('\t')

        return local_linking_number

    linking_number = sum([local_linking_number(vertex) for vertex in crossing_orders])
    return linking_number
