import src.virtual_genus as vg

import sys
import pathlib

source = pathlib.Path(__file__).resolve().parent / sys.argv[1]
target = pathlib.Path(__file__).resolve().parent / sys.argv[2]
assert source.is_file()

with source.open(mode='r') as file:
    with target.open(mode='w') as out:
        for line in file:

            # Deal with comments.
            if line[0] == '#':
                out.write(line)
                break

            fields = line.rstrip(' \n').split(' ')
            properties = fields[1]
            gauss_code = fields[2]

            try:
                g = vg.virtual_genus(gauss_code)
            except ValueError:
                g = None

            if g is None:
                properties += '-'
                genus = 'g-'
            else:
                properties += 'a'
                genus = f'g{g}'
            fields.insert(2, genus)

            fields.append('g' + str(genus))
            line = ' '.join(fields)
            out.write(line + '\n')
