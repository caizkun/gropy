# This script is for testing code.

import sys

from gropy.Gro import Gro


if (len(sys.argv) == 3):
    in_file  = sys.argv[1]
    out_file = sys.argv[2]
else:
    in_file  = 'input.gro'
    out_file = 'output.gro'


test_gro = Gro()
test_gro.read_gro_file(in_file)
test_gro.write_gro_file(out_file)
