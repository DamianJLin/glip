import src.mock_seifert_matrix as msm
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

msms = msm.mock_seifert_matrices(gc, verbose, very_verbose)

print(msms)
