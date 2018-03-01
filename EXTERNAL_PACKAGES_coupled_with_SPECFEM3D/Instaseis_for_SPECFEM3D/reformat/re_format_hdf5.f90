program re_format_hdf5

  use HDF5
  use MPI
  use trim_comments

  implicit none

  ! double precision
  !integer, parameter :: dp = kind(1.d0)
  !integer, parameter :: sp = kind(1.0)
  ! MPI variables
  character(len=500) prname
  integer myrank, nbproc, ierr
  !integer, dimension(MPI_STATUS_SIZE) :: statut

  ! To read from coupling.par file
  character(len=500) meshdirectory
  character(len=500) LOCAL_PATH, TRACT_PATH
  character(len=500) hdf5_file_uncompressed, hdf5_file_compressed
  character(len=500) fields_hdf5_file, meshfem_hdf5_file, hdf5_specfem_output
  character(len=500) fwd_db_path, bwd_db_path, output_dir

  real lat_src, lon_src
  real lat_mesh, lon_mesh, azi_rot
  integer nSpecfem_proc
  character(len=50) frq_min, frq_max, tmin, tmax
  real dt

  ! vr and tr for every time step itime, but for all rec at once
  real, allocatable :: vr(:, :), tr(:, :)
  ! Iterators
  integer irec, itime

  ! Iostatus
  integer ier

  !! HDF5 declarations
  ! Names (file and HDF5 objects)
  ! File name hdf5_file (path to file) declared above (read off par file)
  !character(len=17), parameter :: spec_file = "specfem_dump.hdf5"
  character(len=5),  parameter :: grp_coords = "local"
  character(len=5),  parameter :: grp_fields = "local"
  character(len=19), parameter :: dset_velocity = "local/velocity"
  character(len=12), parameter :: dset_stress = "local/stress"
  character(len=13), parameter :: dset_normals = "local/normals"
  character(len=19), parameter :: dset_spec_d = "displacement"
  character(len=12), parameter :: dset_spec_s = "strain"
  character(len=9),  parameter :: attr_nbrec = "nb_points"
  character(len=12), parameter :: attr_ntime = "nb_timesteps"
  character(len=26), parameter :: attr_nbrec_by_proc = &
                                    "nb_points_per_specfem_proc"
  character(len=23), parameter :: attr_offset_by_proc = &
                                    "offset_per_specfem_proc"

  ! Identifiers
  integer(hid_t) :: coords_file_id         ! File identifier
  integer(hid_t) :: fields_file_id         ! File identifier
  !integer(hid_t) :: spec_file_id           ! File identifier
  integer(hid_t) :: plist_id               ! Property list identifier
  integer(hid_t) :: plist_id2              ! Property list identifier
  !integer(hid_t) :: plist_id3              ! Property list identifier
  integer(hid_t) :: coords_grp_id          ! Group identifier
  integer(hid_t) :: fields_grp_id          ! Group identifier
  !integer(hid_t) :: spec_grp_id            ! Group identifier
  integer(hid_t) :: dset_v_id       ! Dataset identifier
  integer(hid_t) :: dset_s_id       ! Dataset identifier
  !integer(hid_t) :: dset_spec_d_id  ! Dataset identifier
  !integer(hid_t) :: dset_spec_s_id  ! Dataset identifier
  integer(hid_t) :: dset_n_id       ! Dataset identifier
  integer(hid_t) :: dspace_v_id     ! Dataspace identifier
  integer(hid_t) :: dspace_s_id     ! Dataspace identifier
  !integer(hid_t) :: dspace_spec_d_id     ! Dataspace identifier
  !integer(hid_t) :: dspace_spec_s_id     ! Dataspace identifier
  integer(hid_t) :: dspace_n_id     ! Dataspace identifier
  !integer(hid_t) :: dspace_attr_id  ! Dataspace identifier
  integer(hid_t) :: mspace_v_id     ! Memspace identifier
  integer(hid_t) :: mspace_s_id     ! Memspace identifier
  integer(hid_t) :: mspace_n_id     ! Memspace identifier
  integer(hid_t) :: attr_nbrec_id   ! Attribute identifier
  integer(hid_t) :: attr_ntime_id   ! Attribite identifier
  integer(hid_t) :: attr_nbrec_by_proc_id    ! Attribite identifier
  integer(hid_t) :: attr_offset_by_proc_id   ! Attribite identifier

  integer :: error ! Error flag

  real, allocatable :: velocity(:, :, :), stress(:, :, :)     ! Data buffers
  real, allocatable :: normal(:, :)                           ! Data buffers
  integer, allocatable :: nrec_by_proc(:), offset_by_proc(:)  ! Attr buffers
  integer ntime, nbrec                                        ! Attr buffers

  integer :: v_rank = 3, s_rank = 3, n_rank = 2            ! Dataset rank
  integer(hsize_t), dimension(3) :: v_dims, s_dims, n_dims ! Dataset dimensions
  !integer(hsize_t), dimension(3) :: d_dimsf, s_dimsf       ! Dataset dimensions
  integer(hsize_t), dimension(1) :: attr_dims = 1          ! Attribute
                                                           ! dimensions
  integer(hsize_t), dimension(1) :: attr_nbrec_by_proc_dim ! Attribute
                                                           ! dimensions

  integer(hsize_t), dimension(3) :: v_countf, s_countf ! Size of the hyperslab
  integer(hsize_t), dimension(2) :: n_countf           ! in file

  !integer(hsize_t), dimension(2) :: v_countm, s_countm ! Size of the hyperslab
                                                       ! in memory
  integer(hsize_t), dimension(3) :: offsetf            ! Hyperslab offset in f.
  integer(hsize_t), dimension(2) :: n_offsetf          ! Hyperslab offset in f.

  double precision tmp1, tmp2, tmp3     ! To avoid errors

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierr)

  if (myrank == 0) then
    open(10,file='./Inputs_Instaseis_Coupling/coupling.par')
    read(10,'(a)') meshfem_hdf5_file        !! meshfem3D bd points (cartesian)
    read(10,*) nSpecfem_proc                !! number of specfem procs
    read(10,*) lat_src, lon_src             !! axisem source position
    read(10,*) lat_mesh, lon_mesh, azi_rot  !! mesh center position
    read(10,'(a)') meshdirectory            !! mesfem3D results
    read(10,'(a)') LOCAL_PATH               !! specfem database path
    read(10,'(a)') TRACT_PATH               !! path where we save vr&tr
    read(10,'(a)') hdf5_file_uncompressed   !! path of hdf5 file to read
    read(10,'(a)') hdf5_file_compressed     !! path of hdf5 file to read
    read(10,'(a)') hdf5_specfem_output      !! path of hdf5 file (not used here)
    read(10,'(a)') fwd_db_path              !! path of fwd db (not used here)
    read(10,'(a)') bwd_db_path              !! path of bwd db (not used here)
    read(10,'(a)') output_dir               !! path of output dir (not used here)
    read(10,*) frq_min                      !! filter freq (not used here)
    read(10,*) frq_max                      !! filter freq (not used here)
    read(10,*) dt                           !! dt specfem simu (not used here)
    read(10,*) tmin                         !! tmin (not used here)
    read(10,*) tmax                         !! tmax (not used here)
    close(10)

    call trim_string(meshdirectory)
    call trim_string(LOCAL_PATH)
    call trim_string(TRACT_PATH)
    call trim_string(fwd_db_path)
    call trim_string(bwd_db_path)
    fields_hdf5_file = hdf5_file_compressed
    call trim_string(fields_hdf5_file)
    call trim_string(meshfem_hdf5_file)
  endif !! myrank == 0

  call mpi_bcast(nSpecfem_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(LOCAL_PATH, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(TRACT_PATH, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(fields_hdf5_file, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(meshfem_hdf5_file, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  !write(*, *) "!!!!", myrank, meshfem_hdf5_file, len(meshfem_hdf5_file)
  !write(*, *) "!!!!", myrank, fields_hdf5_file, len(fields_hdf5_file)

  allocate(nrec_by_proc(nSpecfem_proc))
  allocate(offset_by_proc(nSpecfem_proc))

  if (myrank == 0) then
    write(*, *) "Opening hdf5 coordinates file..."
  endif

  attr_nbrec_by_proc_dim = nSpecfem_proc
  ! Initialize hdf5 interface
  call h5open_f(error)
  ! Open existing hdf5 file
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  call h5fopen_f(meshfem_hdf5_file, &
                 H5F_ACC_RDONLY_F, coords_file_id, error, &
                 access_prp = plist_id)
  call h5pclose_f(plist_id, error)
  ! Open group to which attributes are attached
  call h5gopen_f(coords_file_id, grp_coords, coords_grp_id, error)
  ! Open attributes
  call h5aopen_f(coords_grp_id, attr_nbrec_by_proc, attr_nbrec_by_proc_id, &
                 error)
  call h5aopen_f(coords_grp_id, attr_offset_by_proc, attr_offset_by_proc_id, &
                 error)

  if (myrank == 0) then
    write(*, *) "Reading hdf5 coordinates file attributes and normals..."
  endif
  ! Read attribute
  call h5aread_f(attr_nbrec_by_proc_id, H5T_NATIVE_INTEGER, &
                 nrec_by_proc, attr_nbrec_by_proc_dim, error)
  call h5aread_f(attr_offset_by_proc_id, H5T_NATIVE_INTEGER, &
                 offset_by_proc, attr_nbrec_by_proc_dim, error)

  n_offsetf = (/0, offset_by_proc(myrank+1)/)
  n_countf = (/3, nrec_by_proc(myrank+1)/)
  n_dims(1) = 3
  n_dims(2) = nrec_by_proc(myrank+1)
  allocate(normal(3, nrec_by_proc(myrank+1)))

  ! Open existing dataset /local/normals
  call h5dopen_f(coords_file_id, dset_normals, dset_n_id, error)

  ! Get datasets' dataspace identifiers
  call h5dget_space_f(dset_n_id, dspace_n_id, error)
  ! Select subset (hyperslab) in the dataset (in the file)
  call h5sselect_hyperslab_f(dspace_n_id, H5S_SELECT_SET_F, &
                             n_offsetf, n_countf, error)
  ! Create memory dataspace
  call h5screate_simple_f(n_rank, n_dims, mspace_n_id, error)
  ! Read data from hyperslab in the file into the hyperslab in memory
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

  call h5dread_f(dset_n_id, H5T_NATIVE_REAL, normal, &
                 n_dims, error, mspace_n_id, dspace_n_id, &
                 xfer_prp = plist_id)


  !! MPC fields output file
  if (myrank == 0) then
    write(*, *) "Opening hdf5 fields output file..."
  endif

  ! Open existing hdf5 file
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id2, error)
  call h5pset_fapl_mpio_f(plist_id2, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  call h5fopen_f(fields_hdf5_file, &
                H5F_ACC_RDONLY_F, fields_file_id, error, &
                access_prp = plist_id2)
  call h5pclose_f(plist_id2, error)

  ! Open group to which attributes are attached
  call h5gopen_f(fields_file_id, grp_fields, fields_grp_id, error)
  ! Open attributes
  call h5aopen_f(fields_grp_id, attr_nbrec, attr_nbrec_id, error)
  call h5aopen_f(fields_grp_id, attr_ntime, attr_ntime_id, error)

  if (myrank == 0) then
    write(*, *) "Reading hdf5 fields output file attributes and datasets..."
  endif

  ! Read attributes
  call h5aread_f(attr_nbrec_id, H5T_NATIVE_INTEGER, &
                 nbrec, attr_dims, error)
  call h5aread_f(attr_ntime_id, H5T_NATIVE_INTEGER, &
                 ntime, attr_dims, error)

  ! Assign dimensions for datasets to read: [icomp, itime, nbrec]
  allocate(velocity(3, ntime, nrec_by_proc(myrank+1)))
  allocate(stress(6, ntime, nrec_by_proc(myrank+1)))

  v_dims(1) = 3
  v_dims(2) = ntime
  v_dims(3) = nrec_by_proc(myrank+1)

  s_dims(1) = 6
  s_dims(2) = ntime
  s_dims(3) = nrec_by_proc(myrank+1)

  v_countf = (/3, ntime, nrec_by_proc(myrank+1)/)
  s_countf = (/6, ntime, nrec_by_proc(myrank+1)/)

  ! Open existing datasets /local/velocity and /local/stress
  call h5dopen_f(fields_file_id, dset_velocity, dset_v_id, error)
  call h5dopen_f(fields_file_id, dset_stress, dset_s_id, error)
  ! Get datasets' dataspace identifiers
  call h5dget_space_f(dset_v_id, dspace_v_id, error)
  call h5dget_space_f(dset_s_id, dspace_s_id, error)

  offsetf  = (/0, 0, offset_by_proc(myrank+1)/)

  ! Select subset (hyperslab) in the dataset (in the file)
  call h5sselect_hyperslab_f(dspace_v_id, H5S_SELECT_SET_F, &
                             offsetf, v_countf, error)
  call h5sselect_hyperslab_f(dspace_s_id, H5S_SELECT_SET_F, &
                             offsetf, s_countf, error)
  ! Create memory dataspace
  call h5screate_simple_f(v_rank, v_dims, mspace_v_id, error)
  call h5screate_simple_f(s_rank, s_dims, mspace_s_id, error)
  ! Read data from hyperslab in the file into the hyperslab in memory
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id2, error)
  call h5pset_dxpl_mpio_f(plist_id2, H5FD_MPIO_INDEPENDENT_F, error)

  call h5dread_f(dset_v_id, H5T_NATIVE_REAL, velocity, &
                 v_dims, error, mspace_v_id, dspace_v_id, &
                 xfer_prp = plist_id2)
  call h5dread_f(dset_s_id, H5T_NATIVE_REAL, stress, &
                 s_dims, error, mspace_s_id, dspace_s_id, &
                 xfer_prp = plist_id2)

  call MPI_Barrier(MPI_COMM_WORLD,ier)


  !! MPC write the specfem binary files
  call create_name_database(prname, myrank, TRACT_PATH)
  write(*, *) "Process", myrank, "opening", TRACT_PATH

  open(111,file=prname(1:len_trim(prname))//'info_NaN',status='unknown', &
     action='write',form='formatted', iostat=ier)

  open(28,file=prname(1:len_trim(prname))//'sol_axisem',status='unknown', &
     action='write',form='unformatted', iostat=ier)

  if (ier /= 0) write(*, *) &
     'error opening', prname(1:len_trim(prname))//'sol_axisem'

  write(*, *) "Process", myrank, &
   "computing tractions and writing out velocities and tractions..."
  allocate(vr(3, nrec_by_proc(myrank+1)), tr(3, nrec_by_proc(myrank+1)))

  do itime = 1, ntime
   do irec = 1, nrec_by_proc(myrank+1)

      vr(1, irec) = velocity(1, itime, irec)
      vr(2, irec) = velocity(2, itime, irec)
      vr(3, irec) = velocity(3, itime, irec)

      tmp1 = dble(stress(1, itime, irec)) * dble(normal(1, irec)) &
         + dble(stress(6, itime, irec)) * dble(normal(2, irec)) &
         + dble(stress(5, itime, irec)) * dble(normal(3, irec))
      tr(1, irec) = real(tmp1)
      tmp2 = dble(stress(6, itime, irec)) * dble(normal(1, irec)) &
         + dble(stress(2, itime, irec)) * dble(normal(2, irec)) &
         + dble(stress(4, itime, irec)) * dble(normal(3, irec))
      tr(2, irec) = real(tmp2)
      tmp3 = dble(stress(5, itime, irec)) * dble(normal(1, irec)) &
         + dble(stress(4, itime, irec)) * dble(normal(2, irec)) &
         + dble(stress(3, itime, irec)) * dble(normal(3, irec))
      tr(3, irec) = real(tmp3)

   enddo !! irec
   !! MPC vr and tr for every time step itime, but for all rec at once
   write(28) vr, tr
  enddo !! itime

  close(111)
  close(28)
  call MPI_Barrier(MPI_COMM_WORLD,ier)

  call h5aclose_f(attr_nbrec_by_proc_id, error)
  call h5aclose_f(attr_offset_by_proc_id, error)
  call h5dclose_f(dset_n_id, error)
  call h5gclose_f(coords_grp_id, error)
  call h5pclose_f(plist_id, error)
  call h5fclose_f(coords_file_id, error)

  ! Close hdf5 groups, datasets, dataspaces attributes...
  call h5aclose_f(attr_nbrec_id, error)
  call h5aclose_f(attr_ntime_id, error)
  call h5sclose_f(dspace_v_id, error)
  call h5sclose_f(dspace_s_id, error)
  call h5sclose_f(mspace_v_id, error)
  call h5sclose_f(mspace_s_id, error)
  call h5pclose_f(plist_id2, error)
  call h5dclose_f(dset_v_id, error)
  call h5dclose_f(dset_s_id, error)
  call h5gclose_f(fields_grp_id, error)
  ! Close the hdf5 file
  call h5fclose_f(fields_file_id, error)

  !if (myrank == 0) then
  !  write(*, *) "Preparing a hdf5 file for specfem dumps..."
  !endif
  !! MPC Create a file that is later going to be filled out by specfem
  ! Setup file access property list with parallel I/O access.
  !call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id3, error)
  !call h5pset_fapl_mpio_f(plist_id3, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! Create the file collectively.
  !call h5fcreate_f(spec_file, H5F_ACC_TRUNC_F, spec_file_id, error, &
  !                 access_prp = plist_id3)
  !call h5pclose_f(plist_id3, error)
  ! Create new group
  !call h5gcreate_f(spec_file_id, grp_coords, spec_grp_id, error)

  !d_dimsf(1) = 3
  !d_dimsf(2) = ntime
  !d_dimsf(3) = nbrec
  !s_dimsf(1) = 6
  !s_dimsf(2) = ntime
  !s_dimsf(3) = nbrec

  ! Create data spaces for the datasets, the dims are the entire arrays
  !call h5screate_simple_f(v_rank, d_dimsf, dspace_spec_d_id, error)
  !call h5screate_simple_f(s_rank, s_dimsf, dspace_spec_s_id, error)
  ! Create datasets with default properties.
  !call h5dcreate_f(spec_grp_id, dset_spec_d, H5T_NATIVE_REAL, &
  !                 dspace_spec_d_id, dset_spec_d_id, error)
  !call h5dcreate_f(spec_grp_id, dset_spec_s, H5T_NATIVE_REAL, &
  !                 dspace_spec_s_id, dset_spec_s_id, error)

  !call h5dclose_f(dset_spec_d_id, error)
  !call h5dclose_f(dset_spec_s_id, error)
  !call h5gclose_f(spec_grp_id, error)
  ! Close the file.
  !call h5fclose_f(spec_file_id, error)

  if (myrank == 0) then
    write(*, *) "Closing the hdf5 interface..."
  endif
  ! Close hdf5 interface
  call h5close_f(error)

  write(*,*) 'Process', myrank, &
    'is done. Number of receivers and timesteps written to binary:', &
    nrec_by_proc(myrank+1), ntime

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  write(*, *) myrank, ierr
  call MPI_FINALIZE(ierr)
  write(*, *) myrank, ierr

  stop

end program re_format_hdf5


!=====================================================================

subroutine create_name_database(prname, iproc, LOCAL_PATH)

! create the name of the database for the mesher and the solver

  implicit none

  integer iproc

! name of the database file
  character(len=500) prname, procname, LOCAL_PATH, clean_LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_')") iproc

! suppress white spaces if any
  clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full name with path
  prname = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

end subroutine create_name_database
