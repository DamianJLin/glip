import sympy as sp
import src.gordon_litherland as gl
import lib.invariant as inv

if __name__ == '__main__':
    import pathlib

    source = pathlib.Path(__file__).resolve().parent / 'table.txt'
    target = pathlib.Path(__file__).resolve().parent / 'out.txt'
    assert source.is_file()

    for g in (0, 1, 2, 3):
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
                        form = gl.gordon_litherland_green(gauss_code)
                    except ValueError:
                        form = 'not alternating'
                    except AssertionError:
                        form = 'asserion error'
                    except:
                        form = 'other error'

                    if fields[2][1] == str(g):
                        out.write(
                            f'{name} {properties} {genus} {gauss_code}\n'
                        )
                        if isinstance(form, str):
                            out.write(
                                f'{form}\n\n'
                            )
                        if not isinstance(form, str):
                            out.write(
                                f'{sp.pretty(form)}, '
                                f'{inv.determinant(form)}, {inv.kobayashi(form)}\n'
                            )
                            out.write('\n')
