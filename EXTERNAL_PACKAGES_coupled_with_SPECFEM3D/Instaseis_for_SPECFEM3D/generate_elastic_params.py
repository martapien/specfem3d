import instaseis
import os
import sys
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

# ------ get path to par files ------
if len(sys.argv) != 2:
    raise ValueError("To run the script need to pass the path to the" \
                     "folder with par files.")
else:
    path_to_par = sys.argv[1]

WORK_DIR = os.getcwd()

# ------ read coupling parameters --------
path_to_coupling_pars = path_to_par + '/coupling.par'
# read from coupling.par
f_cpl = open(path_to_coupling_pars, 'r')
lines_cpl = f_cpl.read().splitlines()
f_cpl.close()
#lines_cpl = list(filter(None, (line.split('#')[0].strip() for line in lines)))

bwd_db_path = lines_cpl[17].rstrip()
coordsfile = os.path.join(WORK_DIR, lines_cpl[1].rstrip())

# ------ read src parameters --------

instaseis.hybrid_get_elastic_params_parallel(coordsfile,
                                             bwd_db_path, source=None)

if rank == 0:
    print("Finished adding parameters to coords file!")
