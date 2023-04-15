import src.mock_seifert_matrix as msm
import lib.invariant as inv
import sys
import sympy as sp

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

if '--symmetric' in sys.argv or '-s' in sys.argv:
    symmetric = True
else:
    symmetric = False

msms = msm.mock_seifert_matrices(
    gauss_code=gc,
    symmetric=symmetric,
    verbose=verbose,
    very_verbose=very_verbose
)


def inv_tuple(m):
    all_invs = []
    all_invs.append(inv.dimension(m))
    all_invs.append(inv.determinant(m))
    all_invs.append(inv.kobayashi(m))
    all_invs.append(inv.mock_alexander(m, for_sort=True))

    return tuple(all_invs)


msms = sorted(msms, key=inv_tuple)

if invariants:
    for m in msms:
        print(f'dim: {sp.pretty(inv.dimension(m))}')
        print(f'det: {sp.pretty(inv.determinant(m))}')
        print(f'kob: {sp.pretty(inv.kobayashi(m))}')
        print(f'alex: {sp.latex(inv.mock_alexander(m).as_expr())}')
        print()
