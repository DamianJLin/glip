import src.gordon_litherland as gl
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

gl.gordon_litherland_green(gc, verbose, very_verbose)
