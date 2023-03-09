import src.gordon_litherland as gl
import lib.invariant as inv
import sys

gc = sys.argv[1]

if '--verbose' in sys.argv or '-v' in sys.argv:
    verbose = True
else:
    verbose = False

if '--very-verbose' in sys.argv or '-vv' in sys.argv:
    very_verbose = True
else:
    very_verbose = False

if '--invariants' in sys.argv or '-i' in sys.argv:
    invariants = True
else:
    invariants = False

form = gl.gordon_litherland(gc, verbose, very_verbose)

if invariants:
    print(f'Det: \t\t{inv.determinant(form)}')
    print(f'Kobayashi: \t{inv.kobayashi(form)}')
