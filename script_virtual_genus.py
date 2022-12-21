import src.virtual_genus as vg

import sys
import pathlib

source = pathlib.Path(__file__).resolve().parent / sys.argv[1]
target = pathlib.Path(__file__).resolve().parent / sys.argv[2]
assert source.is_file()

with source.open(mode='r') as file:
    with target.open(mode='w') as out:
        for line in file:
            fields = line.rstrip(' \n').split(' ')
            gauss_code = line.rstrip(' \n').split(' ')[2]
            genus = vg.virtual_genus(gauss_code)
            fields.append(str(genus))
            line = ' '.join(fields)
            out.write(line + '\n')
