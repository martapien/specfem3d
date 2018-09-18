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
dt = float(lines_cpl[25])
fieldsfile = os.path.join(WORK_DIR, lines_cpl[3].rstrip())
coordsfile = os.path.join(WORK_DIR, lines_cpl[1].rstrip())
writefolder = os.path.join(WORK_DIR, lines_cpl[19].rstrip())
bg_field_file = os.path.join(WORK_DIR, lines_cpl[11].rstrip())


if rank == 0:
    if not os.path.exists(writefolder):
        os.mkdir(writefolder)


# ------ read receivers --------
path_to_receiver_pars = path_to_par + '/receivers.par'
f_rec = open(path_to_receiver_pars, 'r')
receivers_input = f_rec.read().splitlines()

for i in np.arange(len(receivers_input)):
    lat, lon = receivers_input[i].split()
    receiver = instaseis.Receiver(latitude=lat, longitude=lon)
    # ------ generate seismograms ------

    st_hyb = instaseis.hybrid_get_seismograms_parallel(fieldsfile=fieldsfile,
                                                       coordsfile=coordsfile,
                                                       bg_field_file=bg_field_file,
                                                       receiver=receiver,
                                                       bwd_db_path=bwd_db_path, dt=dt)
    # comm.Barrier()

    if rank == 0:
        # ------ write the seismograms to folder ------
        writepath = os.path.join(writefolder, "hybrid_seismogram_%d.mseed" %i)
        st_hyb.write(writepath, format="MSEED")

if rank == 0:
    print("Finished generating hybrid seismograms!")
