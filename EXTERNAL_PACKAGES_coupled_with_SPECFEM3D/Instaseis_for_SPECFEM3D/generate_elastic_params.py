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
lines = f_cpl.readlines()
f_cpl.close()
lines_cpl = list(filter(None, (line.split('#')[0].strip() for line in lines)))

fwd_db_path = lines_cpl[10]
coordsfile = os.path.join(WORK_DIR, lines_cpl[0])

# ------ read src parameters --------
path_to_source_par = path_to_par + '/source.par'
f_src = open(path_to_source_par, 'r')
lines = f_src.readlines()
f_src.close()
lines_src = list(filter(None, (line.split('#')[0].strip() for line in lines)))

src_lat = float(lines_src[0])
src_lon = float(lines_src[1])
depth_in_m = float(lines_src[2])
m_rr = float(lines_src[3])
m_tt = float(lines_src[4])
m_pp = float(lines_src[5])
m_rt = float(lines_src[6])
m_rp = float(lines_src[7])
m_tp = float(lines_src[8])
source = instaseis.Source(
     latitude=src_lat, longitude=src_lon, depth_in_m=depth_in_m,
     m_rr=m_rr, m_tt=m_tt, m_pp=m_pp,
     m_rt=m_rt, m_rp=m_rp, m_tp=m_tp)


instaseis.hybrid_get_elastic_params_parallel(source, coordsfile, 
                                             fwd_db_path)

if rank == 0:
    print("Finished adding parameters to coords file!")
