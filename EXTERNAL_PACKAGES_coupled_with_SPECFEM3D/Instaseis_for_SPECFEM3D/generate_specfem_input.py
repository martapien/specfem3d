#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate input for specfem in HDF5.
"""

import sys
import os
import instaseis
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

# Get working directory
WORK_DIR = os.getcwd()

# get path to par files
if len(sys.argv) != 2:
    raise ValueError("To run the script need to pass the path to the" \
                     "folder with par files.")
else:
    path_to_par = sys.argv[1]

if rank == 0:

    # define paths relative to path_to_par
    path_to_source_par = path_to_par + '/source.par'
    path_to_coupling_pars = path_to_par + '/coupling.par'
    path_to_hdf5 = path_to_par + '/../../gll_coordinates.hdf5'

    # read from coupling.par
    f_cpl = open(path_to_coupling_pars, 'r')
    lines = f_cpl.readlines()
    f_cpl.close()
    lines_cpl = list(filter(None, (line.split('#')[0].strip() for line in lines)))
    cpl_dict = {i: lines_cpl[i] for i in range(0, len(lines_cpl))}

    # Read from source.par file
    f_src = open(path_to_source_par, 'r')
    lines = f_src.readlines()
    f_src.close()
    lines_src = \
        list(filter(None, (line.split('#')[0].strip() for line in lines)))
    src_dict = {i: lines_src[i] for i in range(0, len(lines_src))}
else:
    cpl_dict = None
    src_dict = None
    path_to_hdf5 = None

cpl_dict = comm.bcast(cpl_dict, root=0)
src_dict = comm.bcast(src_dict, root=0)
path_to_hdf5 = comm.bcast(path_to_hdf5, root=0)

outputfile = os.path.join(WORK_DIR, cpl_dict[7])
fwd_db_path = cpl_dict[10]
fmin = cpl_dict[13]
fmax = cpl_dict[14]
dt = float(cpl_dict[15])
tmin = cpl_dict[16]
tmax = cpl_dict[17]

if rank == 0:
    print("Instaseis:")
    print("dt: ", dt)
    print("tmin: ", float(tmin))
    print("tmax: ", float(tmax))
if fmin == "None":
    fmin = None
else:
    fmin = float(fmin)
if fmax == "None":
    fmax = None
else:
    fmax = float(fmax)

if tmin == "None" or tmax == "None":
    time_window = None
else:
    time_window = (float(tmin), float(tmax))


src_lat = float(src_dict[0])
src_lon = float(src_dict[1])
depth_in_m = float(src_dict[2])
m_rr = float(src_dict[3])
m_tt = float(src_dict[4])
m_pp = float(src_dict[5])
m_rt = float(src_dict[6])
m_rp = float(src_dict[7])
m_tp = float(src_dict[8])
source = instaseis.Source(
     latitude=src_lat, longitude=src_lon, depth_in_m=depth_in_m,
     m_rr=m_rr, m_tt=m_tt, m_pp=m_pp,
     m_rt=m_rt, m_rp=m_rp, m_tp=m_tp)

# define dump fields
dumpfields = ("velocity", "stress")
# define dump coordinates
dumpcoords = "local"

# define filter frequencies
if fmin is not None and fmax is not None:
    filter_freqs = (fmin, fmax)
else:
    filter_freqs = None

# the output is generated in local coordinates!
# it later is taken by xreformat!
instaseis.hybrid_generate_output_parallel(path_to_hdf5,
                                          outputfile, fwd_db_path, dt,
                                          source, time_window=time_window,
                                          dumpcoords=dumpcoords,
                                          dumpfields=dumpfields,
                                          remove_source_shift=False)
#                                         filter_freqs=filter_freqs)
comm.Barrier()

if rank == 0:
    print("Finished generate_output_parallel!")
