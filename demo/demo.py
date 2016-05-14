# This script is for demonstration.

import sys

from gropy.Gro import Gro


if (len(sys.argv) == 3):
    in_file  = sys.argv[1]
    out_file = sys.argv[2]
else:
    in_file  = 'input.gro'
    out_file = 'output.gro'


system_gro = Gro()
system_gro.read_gro_file(in_file)
system_gro.sort_residues(['CL', 'SOL'])
system_gro.write_gro_file(out_file)
