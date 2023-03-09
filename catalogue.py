import sys
import functools
import sympy as sp
import src.mock_seifert_matrix as msm
import lib.invariant as inv
import lib.error as er
from progress.bar import IncrementalBar


@functools.total_ordering
class MatDetKob:
    """
    Class containing matrix, determinant and Kobayashi invariant, with
    ordering given by dimension, determinant and Kobayashi invariant heirarchically.
    """

    def __init__(self, mat, det, kob):
        self.mat = mat
        self.det = det
        self.kob = kob

    def __lt__(self, other):
        if sp.shape(self.mat) != sp.shape(other.mat):
            return sp.shape(self.mat) < sp.shape(other.mat)
        elif self.det != other.det:
            return self.det < other.det
        elif self.kob != other.kob:
            return self.kob < other.kob
        elif list(self.mat) != list(other.mat):
            return list(self.mat) < list(other.mat)
        else:
            return False

    def __eq__(self, other):
        return self.det == other.det and self.kob == other.kob and self.max == other.mat


if __name__ == '__main__':
    import pathlib

    symmetric = False
    if '-s' in sys.argv or '--symmetric' in sys.argv:
        symmetric = True

    source = pathlib.Path(__file__).resolve().parent / 'table.txt'
    target = pathlib.Path(__file__).resolve().parent / 'out.txt'
    assert source.is_file()

    genus = (0, 1, 2, 3)
    n_steps = sum(1 for _ in open(source)) * len(genus)
    bar = IncrementalBar('Progress:', max=n_steps, suffix='%(percent)d%%')

    for g in genus:
        with source.open(mode='r') as file:
            with target.open(mode='w' if g == 0 else 'a') as out:

                out.write(f'# ====== Genus {g} ======\n\n')

                for line in file:

                    # Deal with comments.
                    if line[0] == '#':
                        out.write(line)
                        break

                    fields = line.rstrip(' \n').split(' ')
                    name = fields[0]
                    properties = fields[1]
                    genus = fields[2]
                    gauss_code = fields[3]

                    try:
                        m1, m2 = msm.mock_seifert_matrices(gauss_code, symmetric=symmetric)
                    except er.GaussCodeNotAlternatingError:
                        bar.next()
                        continue

                    data_1 = (inv.dimension(m1), inv.determinant(m1),
                              inv.kobayashi(m1), inv.mock_alexander(m1))
                    data_2 = (inv.dimension(m2), inv.determinant(m2),
                              inv.kobayashi(m2), inv.mock_alexander(m2))

                    # Sort.
                    if data_1 < data_2:
                        msms = (m1, m2)
                    else:
                        msms = (m2, m1)

                    if fields[2][1] == str(g):
                        out.write(
                            f'{name} {properties} {genus} {gauss_code}\n'
                        )
                        for m in msms:
                            out.write(
                                f'{sp.pretty(m)},\ndet: {inv.determinant(m)},'
                                f'\nkob: {inv.kobayashi(m)},\n'
                                f'alex: {sp.latex(inv.mock_alexander(m).as_expr())}\n'
                            )
                            out.write('\n')
                    bar.next()
print()
