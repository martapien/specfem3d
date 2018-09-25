!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,iphase, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

  !! MPC for Instaseis-Specfem HDF5 dumps
  use HDF5

  use specfem_par_elastic, only: displ, & !epsilon_trace_over_3, &
          epsilondev_xx, epsilondev_yy, epsilon_trace_new, &
          epsilondev_xy, epsilondev_xz, epsilondev_yz

  use specfem_par, only: it_dsm, it_fk, Veloc_dsm_boundary, Tract_dsm_boundary, &
            Veloc_axisem, Tract_axisem, Tract_axisem_time, myrank, sizeprocs

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE, &
                  INJECTION_TECHNIQUE_TYPE, INSTASEIS_INJECTION_BOX_LOCATION, &
                  RECIPROCITY_AND_KH_INTEGRAL, DT

! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  use specfem_par_coupling, only: npt,nbdglb, NTIME_BETWEEN_FFT, &
     VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, NP_RESAMP, &
     vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

!! comment from Vadim Monteiller, Feb 2017:

! txxbd est calcule dans add_to_compute_stacey_viscoelastic_1.F90 qui lui-meme appelle une subroutine qui se trouve
! dans add_to_compute_stacey_viscoelastic_11.F90. En fait je ne peux pas directement stocker txxbd en memoire
! (sinon il faut 500Go de ram, mais j'ai pas ca). Donc tous les 100 pas de temps je fais une fft
! et je prends la partie du sismo qui m'interesse. C'est a ce moment qu txxbd est rempli. C'est le call suivant qui fait ca:
!
! call store_next_FK_solution( VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f, &
!                    WKS_CMPLX_FOR_FFT, WKS_REAL_FOR_FFT, NF_FOR_STORING, &
!                    NF_FOR_FFT, NTIME_BETWEEN_FFT, NPOW_FOR_FFT, &
!                    vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, npt, it, deltat)

! je prends la solution en frequence : VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f
! et je la sors en temps : vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, pour les 100 pas de temps suivants.

! la subroutine store_next_FK_solution est dans add_to_compute_stacey_viscoelastic_11.F90 ligne 882

! on stocke directement la traction (on fait toujours ca avec DSM et AxiSEM aussi),
! ca evite de stocker 6 composantes de stress, surtout qu'on a des problemes de memoire.


! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
! FK surface
  integer :: ipt,ixglob
  integer :: ii, kk, iim1, iip1, iip2
  double precision :: cs(4), w
  real(kind=CUSTOM_REAL) ::  cs_single(4) !vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
! *********************************************************************************

  !! MPC for Instaseis-Specfem HDF5 dumps
  ! Names (file and HDF5 objects)
  character(len=500) hdf5_file_write
  character(len=500) hdf5_file_read
  character(len=5), parameter :: grp_local = "local"
  character(len=18), parameter :: dset_d = "local/displacement"
  character(len=18), parameter :: dset_s = "local/strain"
  character(len=18), parameter :: dset_gll = "gll_weights"
  character(len=26), parameter :: attr_nbrec_by_proc = &
                                    "nb_points_per_specfem_proc"
  character(len=23), parameter :: attr_offset_by_proc = &
                                    "offset_per_specfem_proc"

  ! Identifiers
  integer(hid_t) :: write_file_id     ! File identifier
  integer(hid_t) :: read_file_id     ! File identifier
  integer(hid_t) :: plist_id          ! Property list identifier
  integer(hid_t) :: read_grp_id      ! Group identifier
  integer(hid_t) :: write_grp_id      ! Group identifier
  integer(hid_t) :: dset_d_id         ! Dataset identifier
  integer(hid_t) :: dset_s_id         ! Dataset identifier
  integer(hid_t) :: dspace_d_id       ! Dataspace identifier
  integer(hid_t) :: dspace_s_id       ! Dataspace identifier
  integer(hid_t) :: mspace_d_id       ! Memspace identifier
  integer(hid_t) :: mspace_s_id       ! Memspace identifier
  integer(hid_t) :: attr_nbrec_by_proc_id    ! Attribite identifier
  integer(hid_t) :: attr_offset_by_proc_id   ! Attribite identifier
  integer(hid_t) :: dset_gll_id         ! Dataset identifier
  integer(hid_t) :: dspace_gll_id       ! Dataspace identifier
  integer(hid_t) :: mspace_gll_id       ! Memspace identifier

  integer :: error ! Error flag

  integer :: d_rank = 2, s_rank = 2, gll_rank = 1 ! Dataset ranks in memory and file

  integer(hsize_t), dimension(2) :: d_dimsm     ! Dataset:
  integer(hsize_t), dimension(2) :: s_dimsm     ! - dimensions in memory
  integer(hsize_t), dimension(1) :: gll_dimsm   ! - dimensions in memory
  integer(hsize_t), dimension(1) :: gll_dimsf   ! - dimensions in file
  integer(hsize_t), dimension(3) :: d_countf    ! Hyperslab size in file
  integer(hsize_t), dimension(3) :: s_countf   ! Hyperslab size in file
  integer(hsize_t), dimension(1) :: gll_countf    ! Hyperslab size in file
  integer(hsize_t), dimension(3) :: d_offsetf   ! Hyperslab offset in f
  integer(hsize_t), dimension(3) :: s_offsetf  ! Hyperslab offset in f
  integer(hsize_t), dimension(1) :: gll_offsetf  ! Hyperslab offset in f
  integer(hsize_t), dimension(1) :: attr_nbrec_by_proc_dim   ! Attribute
                                                             ! dimensions
  integer :: dt_rank = 1                        ! Attribure rank
  integer(hsize_t), dimension(1) :: dt_dims = 1 ! Attribute dimension
  integer(hid_t) :: aspace_dt_id    ! Memspace identifier
  character(len=17), parameter :: attr_dt = "dt"
  integer(hid_t) :: attr_dt_id    ! Attribite identifier

  ! Data and attribute buffers
  real, allocatable :: displ_buf(:, :), strain_buf(:, :), weights_buf(:)
  integer, allocatable :: nrec_by_proc(:), offset_by_proc(:)

  ! And some helpers
  integer :: ipoint, nbrec
  character(len=100) line


  !! CD modif. : begin (implemented by VM) !! For coupling with DSM

  integer :: kaxisem, ip


  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
     ipt = 0
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then

      if (old_DSM_coupling_from_Vadim) then
        if (iphase == 1) then
          if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
            call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
          endif
        endif
      else
        !! MODIFS DE NOBU 2D
      endif

    else if ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) .or. &
            ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
             (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_SOURCE))) then

      if (iphase == 1) then
        call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)

        !! CD CD add this
        if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)
      endif

    else if ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
      if (iphase == 1) then

        !! find indices
        ! example: np_resamp = 2 and it = 1,2,3,4,5,6, ..
        !          --> ii = 1,1,2,2,3,3,..
        ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
        ! example: --> kk = 1,2,1,2,1,2,,..
        kk = it - (ii-1) * NP_RESAMP

        w = dble(kk-1) / dble(NP_RESAMP)

        ! Cubic spline values
        cs(4) = w*w*w/6.d0
        cs(1) = 1.d0/6.d0+w*(w-1.d0)/2.d0-cs(4)
        cs(3) = w+cs(1)-2.d0*cs(4)
        cs(2) = 1.d0-cs(1)-cs(3)-cs(4)

        cs_single(:) = sngl(cs(:))

        iim1 = ii-1
        iip1 = ii+1
        iip2 = ii+2

        do ip=1, npt

           vx_FK(ip) = cs_single(1)* VX_t(ip,iim1) + cs_single(2)* VX_t(ip,ii) + cs_single(3)* VX_t(ip,iip1) + &
                cs_single(4)* VX_t(ip,iip2)
           vy_FK(ip) = cs_single(1)* VY_t(ip,iim1) + cs_single(2)* VY_t(ip,ii) + cs_single(3)* VY_t(ip,iip1) + &
                cs_single(4)* VY_t(ip,iip2)
           vz_FK(ip) = cs_single(1)* VZ_t(ip,iim1) + cs_single(2)* VZ_t(ip,ii) + cs_single(3)* VZ_t(ip,iip1) + &
                cs_single(4)* VZ_t(ip,iip2)
           tx_FK(ip) = cs_single(1)* TX_t(ip,iim1) + cs_single(2)* TX_t(ip,ii) + cs_single(3)* TX_t(ip,iip1) + &
                cs_single(4)* TX_t(ip,iip2)
           ty_FK(ip) = cs_single(1)* TY_t(ip,iim1) + cs_single(2)* TY_t(ip,ii) + cs_single(3)* TY_t(ip,iip1) + &
                cs_single(4)* TY_t(ip,iip2)
           tz_FK(ip) = cs_single(1)* TZ_t(ip,iim1) + cs_single(2)* TZ_t(ip,ii) + cs_single(3)* TZ_t(ip,iip1) + &
                cs_single(4)* TZ_t(ip,iip2)

        enddo
        it_fk=it_fk+1
      endif ! (iphase == 1)
    endif ! (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM)
  endif ! (COUPLE_WITH_INJECTION_TECHNIQUE)

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)

  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. &
     (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
     (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_RECEIVER)) then
     if (myrank == 0) then
       open(10,file='./Inputs_Instaseis_Coupling/coupling.par')
       read(10,'(a)') line
       read(10,'(a)') hdf5_file_read        !! meshfem3D bd points (cartesian)
       read(10,'(a)') line
       read(10,'(a)') hdf5_file_write
       close(10)
     endif !! myrank == 0

     call bcast_all_ch_array(hdf5_file_read,1,500)
     call bcast_all_ch_array(hdf5_file_write,1,500)

    !! MPC read off gll_coordinates hdf5 file the offsets & nbpoints per proc
    allocate(nrec_by_proc(sizeprocs))
    allocate(offset_by_proc(sizeprocs))
    !if (myrank == 0) then
    !  write(*, *) "Opening hdf5 coordinates file and reading attributes..."
    !  write(*,* ) "Time step of specfem simulation:", it
    !endif
    attr_nbrec_by_proc_dim = sizeprocs
    ! Initialize hdf5 interface
    call h5open_f(error)
    ! Open existing hdf5 file
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio(plist_id)

    if (it == 1) then
      call h5fopen_f(trim(hdf5_file_read), &
                     H5F_ACC_RDWR_F, read_file_id, error, &
                     access_prp = plist_id)
    else
      call h5fopen_f(trim(hdf5_file_read), &
                     H5F_ACC_RDONLY_F, read_file_id, error, &
                     access_prp = plist_id)
    endif

    call h5pclose_f(plist_id, error)
    ! Open group to which attributes are attached
    call h5gopen_f(read_file_id, grp_local, read_grp_id, error)
    ! Open attributes
    call h5aopen_f(read_grp_id, attr_nbrec_by_proc, attr_nbrec_by_proc_id, &
                   error)
    call h5aopen_f(read_grp_id, attr_offset_by_proc, attr_offset_by_proc_id, &
                   error)
    ! Read attribute
    call h5aread_f(attr_nbrec_by_proc_id, H5T_NATIVE_INTEGER, &
                   nrec_by_proc, attr_nbrec_by_proc_dim, error)
    call h5aread_f(attr_offset_by_proc_id, H5T_NATIVE_INTEGER, &
                   offset_by_proc, attr_nbrec_by_proc_dim, error)

    allocate(displ_buf(3, nrec_by_proc(myrank+1)))
    allocate(strain_buf(6, nrec_by_proc(myrank+1)))
    allocate(weights_buf(nrec_by_proc(myrank+1)))

  endif

  ipoint = 0
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        vx = veloc(1,iglob)
        vy = veloc(2,iglob)
        vz = veloc(3,iglob)

          !! CD CD !! For coupling with EXTERNAL CODE
          if (COUPLE_WITH_INJECTION_TECHNIQUE) then

            if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then  !! To verify for NOBU version
              vx = vx - Veloc_dsm_boundary(1,it_dsm,igll,iface)
              vy = vy - Veloc_dsm_boundary(2,it_dsm,igll,iface)
              vz = vz - Veloc_dsm_boundary(3,it_dsm,igll,iface)

            else if ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) .or. &
                     ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
                     (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_SOURCE))) then !! VM VM add AxiSEM; MPC add Instaseis
                kaxisem = igll + NGLLSQUARE*(iface - 1)
                vx = vx - Veloc_axisem(1,kaxisem)
                vy = vy - Veloc_axisem(2,kaxisem)
                vz = vz - Veloc_axisem(3,kaxisem)
            endif

          endif

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          if (COUPLE_WITH_INJECTION_TECHNIQUE .and. INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
            ipt = ipt + 1

!! DEBUGVM pour eviter de stocker pour profiler la vitesse de FK
            !vx_FK = vxbd(it_fk,ipt)
            !vy_FK = vybd(it_fk,ipt)
            !vz_FK = vzbd(it_fk,ipt)

            ! sanity check, make sure we are at the right point
            ixglob = nbdglb(ipt)
            !write(*,'(3(i10,1x),3x,  3i3, 3x, i10)') ipt, ixglob, iglob, i,j,k,ispec
            if (iglob /= ixglob) stop 'wrong boundary index for FK coupling'

            vx = vx - vx_FK(ipt)
            vy = vy - vy_FK(ipt)
            vz = vz - vz_FK(ipt)
          endif
! *********************************************************************************

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz

        ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

          !! CD CD !! For coupling with DSM
          if (COUPLE_WITH_INJECTION_TECHNIQUE) then

            if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then    !! To verify for NOBU version
              tx = tx - Tract_dsm_boundary(1,it_dsm,igll,iface)
              ty = ty - Tract_dsm_boundary(2,it_dsm,igll,iface)
              tz = tz - Tract_dsm_boundary(3,it_dsm,igll,iface)

            else if ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) .or. &
                    ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
                     (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_SOURCE))) then
                tx = tx - Tract_axisem(1,kaxisem)
                ty = ty - Tract_axisem(2,kaxisem)
                tz = tz - Tract_axisem(3,kaxisem)
            endif
          endif

          if (COUPLE_WITH_INJECTION_TECHNIQUE .and. &
             (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
             (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_RECEIVER)) then
            !! MPC to dump per proc for coupling with Instaseis
            ipoint = ipoint + 1

            displ_buf(1, ipoint) = displ(1, iglob)
            displ_buf(2, ipoint) = displ(2, iglob)
            displ_buf(3, ipoint) = displ(3, iglob)
            ! debug check for NaN
            !if (displ_buf(1, ipoint) /= displ_buf(1, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "displ_buf(1, ipoint) : ", displ_buf(1, ipoint)
            !endif
            !if (displ_buf(2, ipoint) /= displ_buf(2, ipoint)) then
          !    write(*, *) "Rank: ", myrank, "displ_buf(2, ipoint) : ", displ_buf(2, ipoint)
            !endif
            !if (displ_buf(3, ipoint) /= displ_buf(3, ipoint)) then
          !    write(*, *) "Rank: ", myrank, "displ_buf(3, ipoint) : ", displ_buf(3, ipoint)
            !endif
            strain_buf(1, ipoint) = epsilondev_xx(i,j,k,ispec) &
                + (epsilon_trace_new(i,j,k,ispec) / 3._CUSTOM_REAL)
            strain_buf(2, ipoint) = epsilondev_yy(i,j,k,ispec) &
                + (epsilon_trace_new(i,j,k,ispec) / 3._CUSTOM_REAL)
            strain_buf(3, ipoint) = epsilon_trace_new(i,j,k,ispec) &
                - strain_buf(1, ipoint) &
                - strain_buf(2, ipoint)
            strain_buf(4, ipoint) = epsilondev_yz(i,j,k,ispec)
            strain_buf(5, ipoint) = epsilondev_xz(i,j,k,ispec)
            strain_buf(6, ipoint) = epsilondev_xy(i,j,k,ispec)

            ! debug check for NaN
            !if (strain_buf(1, ipoint) /= strain_buf(1, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(1, ipoint) : ", strain_buf(1, ipoint), ipoint
            !endif
            !if (strain_buf(2, ipoint) /= strain_buf(2, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(2, ipoint) : ", strain_buf(2, ipoint), ipoint
            !endif
            !if (strain_buf(3, ipoint) /= strain_buf(3, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(3, ipoint) : ", strain_buf(3, ipoint), ipoint
            !endif
            !if (strain_buf(4, ipoint) /= strain_buf(4, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(4, ipoint) : ", strain_buf(4, ipoint), ipoint
            !endif
            !if (strain_buf(5, ipoint) /= strain_buf(5, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(5, ipoint) : ", strain_buf(5, ipoint), ipoint
            !endif
            !if (strain_buf(6, ipoint) /= strain_buf(6, ipoint)) then
            !  write(*, *) "Rank: ", myrank, "strain_buf(6, ipoint) : ", strain_buf(6, ipoint), ipoint
            !endif
          endif

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          if (COUPLE_WITH_INJECTION_TECHNIQUE .and. INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
            tx = tx - tx_FK(ipt)
            ty = ty - ty_FK(ipt)
            tz = tz - tz_FK(ipt)
          endif
! *********************************************************************************

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)
        if ((INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
            (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_RECEIVER) .and. &
            (it == 1)) then
              weights_buf(ipoint) = jacobianw
        endif
        ! adds stacey term (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! adjoint simulations
        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          b_absorb_field(1,igll,iface) = tx*jacobianw
          b_absorb_field(2,igll,iface) = ty*jacobianw
          b_absorb_field(3,igll,iface) = tz*jacobianw
        endif !adjoint

          !! CD CD added this
          if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
              write(237) b_absorb_field(1,igll,iface), b_absorb_field(2,igll,iface), b_absorb_field(3,igll,iface)
              write(238) displ(1,iglob), displ(2,iglob), displ(3,iglob)
          endif

      enddo
    endif ! ispec_is_elastic
  enddo

  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. &
     (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) .and. &
     (INSTASEIS_INJECTION_BOX_LOCATION /= INSTASEIS_INJECTION_BOX_LOCATION_RECEIVER)) then

    if (it == 1) then
      !! here we save gll weights!
      gll_offsetf = (/offset_by_proc(myrank+1)/)
      gll_countf = (/nrec_by_proc(myrank+1)/)
      gll_dimsf(1) = sum(nrec_by_proc)
      gll_dimsm(1) = nrec_by_proc(myrank+1)
      !Create data space for the dataset, the dims are the entire arrays
      call h5screate_simple_f(gll_rank, gll_dimsf, dspace_gll_id, error)
      ! Create dataset with default properties.
      call h5dcreate_f(read_grp_id, dset_gll, H5T_NATIVE_REAL, &
                       dspace_gll_id, dset_gll_id, error)
      call h5sselect_hyperslab_f(dspace_gll_id, H5S_SELECT_SET_F, &
                                 gll_offsetf, gll_countf, error)
      ! Create memory dataspace, the dims are the arrays per proc
      call h5screate_simple_f(gll_rank, gll_dimsm, mspace_gll_id, error)
      ! Create property list for independent dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio(plist_id)
      ! Write data
      call h5dwrite_f(dset_gll_id, H5T_NATIVE_REAL, weights_buf, &
          gll_dimsm, error, mspace_gll_id, dspace_gll_id, &
          xfer_prp = plist_id)
      call h5sclose_f(dspace_gll_id, error)
      call h5sclose_f(mspace_gll_id, error)
      call h5dclose_f(dset_gll_id, error)
      call h5pclose_f(plist_id, error)
    endif  !! it==1

    call h5aclose_f(attr_nbrec_by_proc_id, error)
    call h5aclose_f(attr_offset_by_proc_id, error)
    call h5gclose_f(read_grp_id, error)
    call h5fclose_f(read_file_id, error)
    call h5close_f(error)

    !! MPC dump displacement and strain in hdf5
    ! Initialize hdf5 interface
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio(plist_id)
    ! Open the file collectively.
    call h5fopen_f(trim(hdf5_file_write), H5F_ACC_RDWR_F, write_file_id, error, &
                   access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    if (it == 1) then
      call h5gopen_f(write_file_id, grp_local, write_grp_id, error)
      call h5screate_simple_f(dt_rank, dt_dims, aspace_dt_id, error)
      call h5acreate_f(write_grp_id, attr_dt, H5T_NATIVE_DOUBLE, &
                    aspace_dt_id, attr_dt_id, error)
      call h5awrite_f(attr_dt_id, H5T_NATIVE_DOUBLE, DT, &
                      dt_dims, error)
      call h5sclose_f(aspace_dt_id, error)
      call h5aclose_f(attr_dt_id, error)
      call h5gclose_f(write_grp_id, error)

    endif ! if it==1

    ! Open existing datasets
    call h5dopen_f(write_file_id, dset_d, dset_d_id, error)
    call h5dopen_f(write_file_id, dset_s, dset_s_id, error)
    nbrec = sum(nrec_by_proc)

    d_offsetf = (/0, it-1, offset_by_proc(myrank+1)/)
    d_countf = (/3, 1, nrec_by_proc(myrank+1)/)
    d_dimsm(1) = 3
    d_dimsm(2) = nrec_by_proc(myrank+1)

    s_offsetf = (/0, it-1, offset_by_proc(myrank+1)/)
    s_countf = (/6, 1, nrec_by_proc(myrank+1)/)
    s_dimsm(1) = 6
    s_dimsm(2) = nrec_by_proc(myrank+1)

    ! Get datasets' dataspace identifiers
    call h5dget_space_f(dset_d_id, dspace_d_id, error)
    call h5dget_space_f(dset_s_id, dspace_s_id, error)

    call h5sselect_hyperslab_f(dspace_d_id, H5S_SELECT_SET_F, &
                               d_offsetf, d_countf, error)
    call h5sselect_hyperslab_f(dspace_s_id, H5S_SELECT_SET_F, &
                               s_offsetf, s_countf, error)
    ! Create memory dataspace, the dims are the arrays per proc
    call h5screate_simple_f(d_rank, d_dimsm, mspace_d_id, error)
    call h5screate_simple_f(s_rank, s_dimsm, mspace_s_id, error)
    ! Create property list for independent dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio(plist_id)
    ! Write data and attributes
    call h5dwrite_f(dset_d_id, H5T_NATIVE_REAL, displ_buf, &
        d_dimsm, error, mspace_d_id, dspace_d_id, &
        xfer_prp = plist_id)
    call h5dwrite_f(dset_s_id, H5T_NATIVE_REAL, strain_buf, &
        s_dimsm, error, mspace_s_id, dspace_s_id, &
        xfer_prp = plist_id)
    ! Close hdf5 groups, datasets, dataspaces attributes...
    call h5sclose_f(dspace_d_id, error)
    call h5sclose_f(dspace_s_id, error)
    call h5sclose_f(mspace_d_id, error)
    call h5sclose_f(mspace_s_id, error)
    call h5dclose_f(dset_d_id, error)
    call h5dclose_f(dset_s_id, error)
    call h5pclose_f(plist_id, error)
    ! Close the file.
    call h5fclose_f(write_file_id, error)
    ! Close hdf5 interface
    call h5close_f(error)
  endif !! (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
    if (iphase == 1) it_dsm = it_dsm + 1
    !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  endif

  end subroutine compute_stacey_viscoelastic

!
!=====================================================================
!

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_backward(NSPEC_AB, &
                        ibool,iphase, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        ispec_is_elastic,SIMULATION_TYPE, &
                        NSTEP,it,NGLOB_ADJOINT,b_accel, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants
  use specfem_par, only: myrank

  implicit none

  integer :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,NGLOB_ADJOINT
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel

! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  ! reads in absorbing boundary array (when first phase is running)
  ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
  call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        ! adjoint simulations
        b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward

!=============================================================================
!
  !! VM VM & CD CD !! For coupling with external code
!
!=============================================================================

  subroutine read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)

  use constants

  implicit none

  integer igll,it_dsm
  integer iface,num_abs_boundary_faces,i,j
  real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  real(kind=CUSTOM_REAL) :: dsm_boundary_tmp(3,Ntime_step_dsm,NGLLX,NGLLY)

  it_dsm = 1
  do iface = 1,num_abs_boundary_faces

    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll + 1
        read(IIN_veloc_dsm) dsm_boundary_tmp(:,:,i,j)
        Veloc_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
        read(IIN_tract_dsm) dsm_boundary_tmp(:,:,i,j)
        Tract_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
      enddo
    enddo
  enddo

  end subroutine read_dsm_file

!=============================================================================

  subroutine read_axisem_file(Veloc_axisem,Tract_axisem,nb)

  use constants

  implicit none

  integer nb
  real(kind=CUSTOM_REAL) :: Veloc_axisem(3,nb)
  real(kind=CUSTOM_REAL) :: Tract_axisem(3,nb)

  read(IIN_veloc_dsm) Veloc_axisem, Tract_axisem

  end subroutine read_axisem_file

!
!=====================================================================
!

! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field,Mesh_pointer)

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  implicit none

! communication overlap
  integer :: iphase

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

  !! For coupling with DSM

!! DK DK beware:: these two arrays are automatic arrays, not subroutine arguments, thus allocated every time
!! DK DK beware:: this routine is called, and erased when the routine exits; that is very strange
! real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
! real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  integer :: it_dsm

  it_dsm = 1 !! initialize dsm iterator
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    if (old_DSM_coupling_from_Vadim) then
      if (iphase == 1) then
        if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
          stop 'DK DK old_DSM_coupling_from_Vadim support is discontinued'
!         call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
        endif
      endif
    else
      !! MODIFS DE NOBU 2D
    endif
  endif

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array (when first phase is running)
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
  endif !adjoint

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_field)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  !! CD CD : begin
  !! For coupling with DSM
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
    if (iphase == 1) then
      it_dsm = it_dsm + 1
    endif
  endif
  !! CD CD : end

  end subroutine compute_stacey_viscoelastic_GPU


