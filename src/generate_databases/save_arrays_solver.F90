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

! for external mesh

  subroutine save_arrays_solver_ext_mesh(nspec,nglob,APPROXIMATE_OCEAN_LOAD,ibool, &
                    num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    SAVE_MESH_FILES,ANISOTROPY)

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,NGLLSQUARE,IMAIN,IOUT, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    SIMULATION_TYPE,SAVE_FORWARD,mask_ibool_interior_domain, &
    STACEY_ABSORBING_CONDITIONS,USE_MESH_COLORING_GPU

  ! PML
  use generate_databases_par, only: PML_CONDITIONS, &
    nspec_cpml,CPML_width_x,CPML_width_y,CPML_width_z,CPML_to_spec, &
    CPML_regions,is_CPML,min_distance_between_CPML_parameter,nspec_cpml_tot, &
    d_store_x,d_store_y,d_store_z,k_store_x,k_store_y,k_store_z, &
    alpha_store_x,alpha_store_y,alpha_store_z, &
    nglob_interface_PML_acoustic,points_interface_PML_acoustic, &
    nglob_interface_PML_elastic,points_interface_PML_elastic

  ! mesh surface
  use generate_databases_par, only: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh, &
    nfaces_surface

  use create_regions_mesh_ext_par

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  integer :: nspec,nglob
  ! ocean load
  logical :: APPROXIMATE_OCEAN_LOAD
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  ! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  logical :: SAVE_MESH_FILES
  logical :: ANISOTROPY

  ! local parameters
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  integer :: ier,i
  character(len=MAX_STRING_LEN) :: filename

  ! saves mesh file proc***_external_mesh.bin
  filename = prname(1:len_trim(prname))//'external_mesh.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening database proc######_external_mesh.bin'

  write(IOUT) nspec
  write(IOUT) nglob

  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

  write(IOUT) xixstore
  write(IOUT) xiystore
  write(IOUT) xizstore
  write(IOUT) etaxstore
  write(IOUT) etaystore
  write(IOUT) etazstore
  write(IOUT) gammaxstore
  write(IOUT) gammaystore
  write(IOUT) gammazstore
  write(IOUT) jacobianstore

  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) ispec_is_acoustic
  write(IOUT) ispec_is_elastic
  write(IOUT) ispec_is_poroelastic

! acoustic
  if (ACOUSTIC_SIMULATION) then
    write(IOUT) rmass_acoustic
  endif

! this array is needed for acoustic simulations but also for elastic simulations with CPML,
! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  write(IOUT) rhostore

! elastic
  if (ELASTIC_SIMULATION) then
    write(IOUT) rmass
    if (APPROXIMATE_OCEAN_LOAD) then
      write(IOUT) rmass_ocean_load
    endif
    !pll Stacey
    write(IOUT) rho_vp
    write(IOUT) rho_vs
  endif

! poroelastic
  if (POROELASTIC_SIMULATION) then
    write(IOUT) rmass_solid_poroelastic
    write(IOUT) rmass_fluid_poroelastic
    write(IOUT) rhoarraystore
    write(IOUT) kappaarraystore
    write(IOUT) etastore
    write(IOUT) tortstore
    write(IOUT) permstore
    write(IOUT) phistore
    write(IOUT) rho_vpI
    write(IOUT) rho_vpII
    write(IOUT) rho_vsI
  endif

! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    write(IOUT) nspec_cpml
    write(IOUT) CPML_width_x
    write(IOUT) CPML_width_y
    write(IOUT) CPML_width_z
    write(IOUT) min_distance_between_CPML_parameter
    if (nspec_cpml > 0) then
      write(IOUT) CPML_regions
      write(IOUT) CPML_to_spec
      write(IOUT) is_CPML
      write(IOUT) d_store_x
      write(IOUT) d_store_y
      write(IOUT) d_store_z
      write(IOUT) k_store_x
      write(IOUT) k_store_y
      write(IOUT) k_store_z
      write(IOUT) alpha_store_x
      write(IOUT) alpha_store_y
      write(IOUT) alpha_store_z
      ! --------------------------------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior computational domain
      ! --------------------------------------------------------------------------------------------
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        write(IOUT) nglob_interface_PML_acoustic
        write(IOUT) nglob_interface_PML_elastic
        if (nglob_interface_PML_acoustic > 0) write(IOUT) points_interface_PML_acoustic
        if (nglob_interface_PML_elastic > 0)  write(IOUT) points_interface_PML_elastic
      endif
    endif
  endif

! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  if (num_abs_boundary_faces > 0) then
    write(IOUT) abs_boundary_ispec
    write(IOUT) abs_boundary_ijk
    write(IOUT) abs_boundary_jacobian2Dw
    write(IOUT) abs_boundary_normal
    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        write(IOUT) rmassx
        write(IOUT) rmassy
        write(IOUT) rmassz
      endif
      if (ACOUSTIC_SIMULATION) then
        write(IOUT) rmassz_acoustic
      endif
    endif
  endif

  write(IOUT) nspec2D_xmin
  write(IOUT) nspec2D_xmax
  write(IOUT) nspec2D_ymin
  write(IOUT) nspec2D_ymax
  write(IOUT) NSPEC2D_BOTTOM
  write(IOUT) NSPEC2D_TOP
  write(IOUT) ibelm_xmin
  write(IOUT) ibelm_xmax
  write(IOUT) ibelm_ymin
  write(IOUT) ibelm_ymax
  write(IOUT) ibelm_bottom
  write(IOUT) ibelm_top

! free surface
  write(IOUT) num_free_surface_faces
  if (num_free_surface_faces > 0) then
    write(IOUT) free_surface_ispec
    write(IOUT) free_surface_ijk
    write(IOUT) free_surface_jacobian2Dw
    write(IOUT) free_surface_normal
  endif

! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  if (num_coupling_ac_el_faces > 0) then
    write(IOUT) coupling_ac_el_ispec
    write(IOUT) coupling_ac_el_ijk
    write(IOUT) coupling_ac_el_jacobian2Dw
    write(IOUT) coupling_ac_el_normal
  endif

! acoustic-poroelastic coupling surface
  write(IOUT) num_coupling_ac_po_faces
  if (num_coupling_ac_po_faces > 0) then
    write(IOUT) coupling_ac_po_ispec
    write(IOUT) coupling_ac_po_ijk
    write(IOUT) coupling_ac_po_jacobian2Dw
    write(IOUT) coupling_ac_po_normal
  endif

! elastic-poroelastic coupling surface
  write(IOUT) num_coupling_el_po_faces
  if (num_coupling_el_po_faces > 0) then
    write(IOUT) coupling_el_po_ispec
    write(IOUT) coupling_po_el_ispec
    write(IOUT) coupling_el_po_ijk
    write(IOUT) coupling_po_el_ijk
    write(IOUT) coupling_el_po_jacobian2Dw
    write(IOUT) coupling_el_po_normal
  endif

  !MPI interfaces
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))

  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo

  write(IOUT) num_interfaces_ext_mesh
  if (num_interfaces_ext_mesh > 0) then
    write(IOUT) max_nibool_interfaces_ext_mesh
    write(IOUT) my_neighbors_ext_mesh
    write(IOUT) nibool_interfaces_ext_mesh
    write(IOUT) ibool_interfaces_ext_mesh_dummy
  endif

! anisotropy
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    write(IOUT) c11store
    write(IOUT) c12store
    write(IOUT) c13store
    write(IOUT) c14store
    write(IOUT) c15store
    write(IOUT) c16store
    write(IOUT) c22store
    write(IOUT) c23store
    write(IOUT) c24store
    write(IOUT) c25store
    write(IOUT) c26store
    write(IOUT) c33store
    write(IOUT) c34store
    write(IOUT) c35store
    write(IOUT) c36store
    write(IOUT) c44store
    write(IOUT) c45store
    write(IOUT) c46store
    write(IOUT) c55store
    write(IOUT) c56store
    write(IOUT) c66store
  endif

! inner/outer elements
  write(IOUT) ispec_is_inner

  if (ACOUSTIC_SIMULATION) then
    write(IOUT) nspec_inner_acoustic,nspec_outer_acoustic
    write(IOUT) num_phase_ispec_acoustic
    if (num_phase_ispec_acoustic > 0) write(IOUT) phase_ispec_inner_acoustic
  endif

  if (ELASTIC_SIMULATION) then
    write(IOUT) nspec_inner_elastic,nspec_outer_elastic
    write(IOUT) num_phase_ispec_elastic
    if (num_phase_ispec_elastic > 0) write(IOUT) phase_ispec_inner_elastic
  endif

  if (POROELASTIC_SIMULATION) then
    write(IOUT) nspec_inner_poroelastic,nspec_outer_poroelastic
    write(IOUT) num_phase_ispec_poroelastic
    if (num_phase_ispec_poroelastic > 0) write(IOUT) phase_ispec_inner_poroelastic
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      write(IOUT) num_colors_outer_acoustic,num_colors_inner_acoustic
      write(IOUT) num_elem_colors_acoustic
    endif
    if (ELASTIC_SIMULATION) then
      write(IOUT) num_colors_outer_elastic,num_colors_inner_elastic
      write(IOUT) num_elem_colors_elastic
    endif
  endif

  ! surface points
  write(IOUT) nfaces_surface
  write(IOUT) ispec_is_surface_external_mesh
  write(IOUT) iglob_is_surface_external_mesh

  close(IOUT)

  ! stores arrays in binary files
  if (SAVE_MESH_FILES) call save_arrays_solver_files(nspec,nglob,ibool)

  ! if SAVE_MESH_FILES is true then the files have already been saved, no need to save them again
  if ((COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) .and. .not. SAVE_MESH_FILES) then
    call save_arrays_solver_files(nspec,nglob,ibool)
  endif

  ! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier)
  if (ier /= 0) stop 'error deallocating array ibool_interfaces_ext_mesh_dummy'

  if (nspec_cpml_tot > 0) then
     deallocate(CPML_to_spec,stat=ier); if (ier /= 0) stop 'error deallocating array CPML_to_spec'
     deallocate(CPML_regions,stat=ier); if (ier /= 0) stop 'error deallocating array CPML_regions'
     deallocate(is_CPML,stat=ier); if (ier /= 0) stop 'error deallocating array is_CPML'
  endif

  if (PML_CONDITIONS) then
     deallocate(d_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_x'
     deallocate(d_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_y'
     deallocate(d_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_z'
     deallocate(k_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_x'
     deallocate(k_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_y'
     deallocate(k_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_z'
     deallocate(alpha_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_x'
     deallocate(alpha_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_y'
     deallocate(alpha_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_z'
     if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
       deallocate(mask_ibool_interior_domain,stat=ier)
       if (ier /= 0) stop 'error deallocating array mask_ibool_interior_domain'

       if (nglob_interface_PML_acoustic > 0) then
         deallocate(points_interface_PML_acoustic,stat=ier)
         if (ier /= 0) stop 'error deallocating array points_interface_PML_acoustic'
       endif

       if (nglob_interface_PML_elastic > 0) then
         deallocate(points_interface_PML_elastic,stat=ier)
         if (ier /= 0) stop 'error deallocating array points_interface_PML_elastic'
       endif
     endif
  endif

  end subroutine save_arrays_solver_ext_mesh

!
!-------------------------------------------------------------------------------------------------
!
  subroutine save_arrays_solver_files(nspec,nglob,ibool)

 !! MPC for instaseis input
  use HDF5

  use generate_databases_par, only: myrank, sizeprocs, &
        NGLLX,NGLLY,NGLLZ,NGLLSQUARE,IMAIN,IOUT,FOUR_THIRDS

  ! MPI interfaces
  use generate_databases_par, only: nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,num_interfaces_ext_mesh


  use create_regions_mesh_ext_par

  use constants, only: INJECTION_TECHNIQUE_IS_INSTASEIS
  use shared_parameters, only: NPROC, COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH, &
    INJECTION_TECHNIQUE_TYPE, NSTEP

  implicit none

  integer :: nspec,nglob
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
  integer,dimension(:),allocatable :: v_tmp_i
  integer :: ier,i,j,k
  integer, dimension(:), allocatable :: iglob_tmp
  integer :: iface, igll, ispec, iglob, inum, num_points
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  character(len=MAX_STRING_LEN) :: filename

  !----------------------------------------------------------------------
  ! mostly for free-surface and coupling surfaces
  logical,parameter :: SAVE_MESH_FILES_ADDITIONAL = .true.


  !! MPC HDF5 declarations for Instaseis-Specfem HDF5 dumps
  ! Names (file and HDF5 objects)
  character(len=20), parameter :: hdf5_file = "gll_coordinates.hdf5" ! File name
  character(len=5),  parameter :: grp_local = "local"
  character(len=18), parameter :: grp_params = "elastic_parameters"
  character(len=18), parameter :: dset_coords = "coordinates"
  character(len=18), parameter :: dset_normals = "normals"
  character(len=24), parameter :: dset_kappa = "kappa"
  character(len=21), parameter :: dset_mu = "mu"
  character(len=22), parameter :: dset_rho = "rho"
  character(len=22), parameter :: attr_nbrec = "nb_points"
  character(len=26), parameter :: attr_spec = "nb_points_per_specfem_proc"
  character(len=23), parameter :: attr_spec2 = "offset_per_specfem_proc"

  ! Identifiers
  integer(hid_t) :: file_id           ! File identifier
  integer(hid_t) :: plist_id          ! Property list identifier
  integer(hid_t) :: grp_local_id      ! Group identifier
  integer(hid_t) :: grp_params_id     ! Group identifier
  integer(hid_t) :: dset_coords_id    ! Dataset identifier
  integer(hid_t) :: dset_normals_id   ! Dataset identifier
  integer(hid_t) :: dset_kappa_id     ! Dataset identifier
  integer(hid_t) :: dset_mu_id        ! Dataset identifier
  integer(hid_t) :: dset_rho_id       ! Dataset identifier
  integer(hid_t) :: dspace_coords_id  ! Dataspace identifier
  integer(hid_t) :: dspace_normals_id ! Dataspace identifier
  integer(hid_t) :: dspace_kappa_id   ! Dataspace identifier
  integer(hid_t) :: dspace_mu_id      ! Dataspace identifier
  integer(hid_t) :: dspace_rho_id     ! Dataspace identifier
  integer(hid_t) :: mspace_coords_id  ! Memspace identifier
  integer(hid_t) :: mspace_normals_id ! Memspace identifier
  integer(hid_t) :: mspace_kappa_id   ! Memspace identifier
  integer(hid_t) :: mspace_mu_id      ! Memspace identifier
  integer(hid_t) :: mspace_rho_id     ! Memspace identifier
  integer(hid_t) :: attr_nbrec_id     ! Attribute identifier
  integer(hid_t) :: attr_spec_id      ! Attribute identifier
  integer(hid_t) :: attr_spec2_id      ! Attribute identifier
  integer(hid_t) :: aspace_nbrec_id   ! Attribute dataspace identifier
  integer(hid_t) :: aspace_spec_id    ! Attribute dataspace identifier
  integer(hid_t) :: aspace_spec2_id    ! Attribute dataspace identifier

  integer :: error ! Error flag

  integer :: coords_rank = 2, params_rank = 1       ! Dataset ranks
                                                    ! in memory and file
  integer :: nbrec_rank = 1, spec_rank = 1          ! Attribure rank
  integer(hsize_t), dimension(2) :: coords_dimsf    ! Dataset
  integer(hsize_t), dimension(2) :: params_dimsf    ! dimensions in file
  integer(hsize_t), dimension(2) :: coords_dimsm    ! Dataset
  integer(hsize_t), dimension(2) :: params_dimsm    ! dimensions in memory
  integer(hsize_t), dimension(1) :: nbrec_dims      ! Attribute dimension
  integer(hsize_t), dimension(1) :: spec_dims      ! Attribute dimension
  integer(hsize_t), dimension(2) :: coords_countf   ! Hyperslab size in file
  integer(hsize_t), dimension(1) :: params_countf   ! Hyperslab size in file
  integer(hsize_t), dimension(2) :: coords_offsetf  ! Hyperslab offset in f
  integer(hsize_t), dimension(1) :: params_offsetf  ! Hyperslab offset in f

  ! Data and attribute buffers
  real, allocatable :: coords_buf(:, :), normals_buf(:, :), &
                       kappa_buf(:), mu_buf(:), rho_buf(:)
  integer :: nbrec, ntime

  ! And some helpers
  integer, allocatable :: offset_per_proc(:), nb_gll_per_proc(:)
  integer :: nb_gll_myrank, ipoint

  !! HDF5 declarations for specfem_dump file
  ! Names (file and HDF5 objects)
  character(len=17), parameter :: spec_file = "specfem_dump.hdf5"
  character(len=5),  parameter :: grp_coords = "local"
  character(len=19), parameter :: dset_spec_d = "displacement"
  character(len=12), parameter :: dset_spec_s = "strain"

  ! Identifiers
  integer(hid_t) :: spec_file_id           ! File identifier
  integer(hid_t) :: spec_grp_id            ! Group identifier
  integer(hid_t) :: dset_spec_d_id  ! Dataset identifier
  integer(hid_t) :: dset_spec_s_id  ! Dataset identifier
  integer(hid_t) :: dspace_spec_d_id     ! Dataspace identifier
  integer(hid_t) :: dspace_spec_s_id     ! Dataspace identifier

  integer(hsize_t), dimension(3) :: d_dimsf, s_dimsf       ! Dataset dimensions
  integer :: v_rank = 3, s_rank = 3           ! Dataset rank


  if (myrank == 0) then
    write(IMAIN,*) '     saving mesh files for AVS, OpenDX, Paraview'
    call flush_IMAIN()
  endif

  ! mesh arrays used for example in combine_vol_data.f90
  !--- x coordinate
  open(unit=IOUT,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file x.bin'
  write(IOUT) xstore_dummy
  close(IOUT)

  !--- y coordinate
  open(unit=IOUT,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file y.bin'
  write(IOUT) ystore_dummy
  close(IOUT)

  !--- z coordinate
  open(unit=IOUT,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file z.bin'
  write(IOUT) zstore_dummy
  close(IOUT)

  ! ibool
  open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file ibool.bin'
  write(IOUT) ibool
  close(IOUT)

  allocate( v_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if (ier /= 0) stop 'error allocating array '

  ! vp (for checking the mesh and model)
  !minimum = minval( abs(rho_vp) )
  !if (minimum(1) /= 0.0) then
  !  v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
  !else
  !  v_tmp = 0.0
  !endif
  v_tmp = 0.0
  where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file vp.bin'
  write(IOUT) v_tmp
  close(IOUT)

  ! vp values - VTK file output
  filename = prname(1:len_trim(prname))//'vp'
  call write_VTK_data_gll_cr(nspec,nglob, &
                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                      v_tmp,filename)


  ! vs (for checking the mesh and model)
  !minimum = minval( abs(rho_vs) )
  !if (minimum(1) /= 0.0) then
  !  v_tmp = mustore / rho_vs
  !else
  !  v_tmp = 0.0
  !endif
  v_tmp = 0.0
  where( rho_vs /= 0._CUSTOM_REAL )  v_tmp = mustore / rho_vs
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file vs.bin'
  write(IOUT) v_tmp
  close(IOUT)

  ! vs values - VTK file output
  filename = prname(1:len_trim(prname))//'vs'
  call write_VTK_data_gll_cr(nspec,nglob, &
                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                      v_tmp,filename)

  ! outputs density model for check
  v_tmp = 0.0
  where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)
  open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file rho.bin'
  write(IOUT) v_tmp
  close(IOUT)

  ! attenuation
  ! shear attenuation Qmu
  open(unit=IOUT,file=prname(1:len_trim(prname))//'qmu.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file qmu.bin'
  write(IOUT) qmu_attenuation_store
  close(IOUT)

  ! shear attenuation - VTK file output
  filename = prname(1:len_trim(prname))//'qmu'
  call write_VTK_data_gll_cr(nspec,nglob, &
                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                      qmu_attenuation_store,filename)

  ! bulk attenuation Qkappa
  open(unit=IOUT,file=prname(1:len_trim(prname))//'qkappa.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file qkappa.bin'
  write(IOUT) qkappa_attenuation_store
  close(IOUT)

  ! bulk attenuation - VTK file output
  filename = prname(1:len_trim(prname))//'qkappa'
  call write_VTK_data_gll_cr(nspec,nglob, &
                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                      qkappa_attenuation_store,filename)

  ! frees temporary array
  deallocate(v_tmp)

  ! additional VTK file output
  if (SAVE_MESH_FILES_ADDITIONAL) then
    ! user output
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '     saving additonal mesh files with surface/coupling points'
      call flush_IMAIN()
    endif

    ! saves free surface points
    if (num_free_surface_faces > 0) then
      ! saves free surface interface points
      allocate( iglob_tmp(NGLLSQUARE*num_free_surface_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_free_surface_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(free_surface_ijk(1,j,i), &
                                  free_surface_ijk(2,j,i), &
                                  free_surface_ijk(3,j,i), &
                                  free_surface_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'free_surface'
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_free_surface_faces, &
                        filename)

      deallocate(iglob_tmp)
    endif

    ! acoustic-elastic domains
    if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
      ! saves points on acoustic-elastic coupling interface
      num_points = NGLLSQUARE*num_coupling_ac_el_faces
      allocate( iglob_tmp(num_points),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i = 1,num_coupling_ac_el_faces
        do j = 1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_ac_el_ijk(1,j,i), &
                                  coupling_ac_el_ijk(2,j,i), &
                                  coupling_ac_el_ijk(3,j,i), &
                                  coupling_ac_el_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_acoustic_elastic'
      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 iglob_tmp,num_points,filename)

      ! saves acoustic/elastic flag
      allocate(v_tmp_i(nspec),stat=ier)
      if (ier /= 0) stop 'error allocating array v_tmp_i'
      do i = 1,nspec
        if (ispec_is_acoustic(i)) then
          v_tmp_i(i) = 1
        else if (ispec_is_elastic(i)) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'acoustic_elastic_flag'
      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                 v_tmp_i,filename)

      deallocate(iglob_tmp,v_tmp_i)
    endif !if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION )

    ! acoustic-poroelastic domains
    if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION) then
      ! saves points on acoustic-poroelastic coupling interface
      num_points = NGLLSQUARE*num_coupling_ac_po_faces
      allocate( iglob_tmp(num_points),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i = 1,num_coupling_ac_po_faces
        do j = 1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_ac_po_ijk(1,j,i), &
                                  coupling_ac_po_ijk(2,j,i), &
                                  coupling_ac_po_ijk(3,j,i), &
                                  coupling_ac_po_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_acoustic_poroelastic'
      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 iglob_tmp,num_points,filename)

      ! saves acoustic/poroelastic flag
      allocate(v_tmp_i(nspec),stat=ier)
      if (ier /= 0) stop 'error allocating array v_tmp_i'
      do i = 1,nspec
        if (ispec_is_acoustic(i)) then
          v_tmp_i(i) = 1
        else if (ispec_is_poroelastic(i)) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'acoustic_poroelastic_flag'
      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                 v_tmp_i,filename)

      deallocate(v_tmp_i,iglob_tmp)
    endif !if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION )

    ! elastic-poroelastic domains
    if (ELASTIC_SIMULATION .and. POROELASTIC_SIMULATION) then
      ! saves points on elastic-poroelastic coupling interface
      num_points = NGLLSQUARE*num_coupling_el_po_faces
      allocate( iglob_tmp(num_points),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i = 1,num_coupling_el_po_faces
        do j = 1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_el_po_ijk(1,j,i), &
                                  coupling_el_po_ijk(2,j,i), &
                                  coupling_el_po_ijk(3,j,i), &
                                  coupling_el_po_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_elastic_poroelastic'
      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 iglob_tmp,num_points,filename)

      ! saves elastic/poroelastic flag
      allocate(v_tmp_i(nspec),stat=ier)
      if (ier /= 0) stop 'error allocating array v_tmp_i'
      do i=1,nspec
        if (ispec_is_elastic(i)) then
          v_tmp_i(i) = 1
        else if (ispec_is_poroelastic(i)) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'elastic_poroelastic_flag'
      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                 v_tmp_i,filename)

      deallocate(v_tmp_i,iglob_tmp)
    endif !if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION

    ! MPI
    if (NPROC > 1) then
      ! saves MPI interface points
      num_points = sum(nibool_interfaces_ext_mesh(1:num_interfaces_ext_mesh))
      allocate( iglob_tmp(num_points),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i = 1,num_interfaces_ext_mesh
        do j = 1, nibool_interfaces_ext_mesh(i)
          inum = inum + 1
          iglob_tmp(inum) = ibool_interfaces_ext_mesh(j,i)
        enddo
      enddo

      filename = prname(1:len_trim(prname))//'MPI_points'
      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 iglob_tmp,num_points,filename)
      deallocate(iglob_tmp)
    endif ! NPROC > 1
  endif  !if (SAVE_MESH_FILES_ADDITIONAL)

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
    !if (num_abs_boundary_faces > 0) then
    filename = prname(1:len_trim(prname))//'absorb_dsm'
    open(IOUT,file=filename(1:len_trim(filename)),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file absorb_dsm'
    write(IOUT) num_abs_boundary_faces
    write(IOUT) abs_boundary_ispec
    write(IOUT) abs_boundary_ijk
    write(IOUT) abs_boundary_jacobian2Dw
    write(IOUT) abs_boundary_normal
    close(IOUT)

    filename = prname(1:len_trim(prname))//'inner'
    open(IOUT,file=filename(1:len_trim(filename)),status='unknown',form='unformatted',iostat=ier)
    write(IOUT) ispec_is_inner
    write(IOUT) ispec_is_elastic
    close(IOUT)

    !! MPC GLL points in local coordinates,kappa, mu and rho

    nb_gll_myrank = 0
    do iface = 1,num_abs_boundary_faces
       ispec = abs_boundary_ispec(iface)
       if (ispec_is_elastic(ispec)) then
          do igll = 1,NGLLSQUARE
             nb_gll_myrank = nb_gll_myrank + 1
          enddo
       endif
    enddo

    allocate(nb_gll_per_proc(sizeprocs))
    allocate(offset_per_proc(sizeprocs))

    call synchronize_all()

    if (myrank > 0) then
      call send_i_t(nb_gll_myrank,1,0)
    else
      nb_gll_per_proc(1) = nb_gll_myrank
      do i=2, sizeprocs
        call recv_i_t(nb_gll_per_proc(i),1,i-1)
      enddo

      offset_per_proc(1) = 0
      do i=2, sizeprocs
        offset_per_proc(i) = sum(nb_gll_per_proc(1:i-1))
      enddo
    endif  !! myrank > 0
    call bcast_all_i(nb_gll_per_proc, sizeprocs)
    call bcast_all_i(offset_per_proc, sizeprocs)


    call synchronize_all()

    allocate(coords_buf(3, nb_gll_per_proc(myrank+1)))
    allocate(normals_buf(3, nb_gll_per_proc(myrank+1)))
    allocate(kappa_buf(nb_gll_per_proc(myrank+1)))
    allocate(mu_buf(nb_gll_per_proc(myrank+1)))
    allocate(rho_buf(nb_gll_per_proc(myrank+1)))

    !! VM VM write an ascii file for instaseis input
    filename = prname(1:len_trim(prname))//'normal.txt'
    open(IOUT,file=filename(1:len_trim(filename)),status='unknown',iostat=ier)
    write(IOUT, *) ' number of points :', num_abs_boundary_faces*NGLLSQUARE

    ipoint = 0
    do iface = 1,num_abs_boundary_faces
       ispec = abs_boundary_ispec(iface)
       if (ispec_is_elastic(ispec)) then
          do igll = 1,NGLLSQUARE

             ! gets local indices for GLL point
             i = abs_boundary_ijk(1,igll,iface)
             j = abs_boundary_ijk(2,igll,iface)
             k = abs_boundary_ijk(3,igll,iface)

             iglob = ibool(i,j,k,ispec)
             ipoint = ipoint + 1

             coords_buf(1, ipoint) = xstore_dummy(iglob)
             coords_buf(2, ipoint) = ystore_dummy(iglob)
             coords_buf(3, ipoint) = zstore_dummy(iglob)

             normals_buf(1, ipoint) = abs_boundary_normal(1,igll,iface)
             normals_buf(2, ipoint) = abs_boundary_normal(2,igll,iface)
             normals_buf(3, ipoint) = abs_boundary_normal(3,igll,iface)

             write(IOUT,'(6f25.10)') xstore_dummy(iglob), ystore_dummy(iglob), &
                zstore_dummy(iglob), normals_buf(1, ipoint), &
                normals_buf(2, ipoint), normals_buf(3, ipoint)
             kappa_buf(ipoint) = kappastore(i,j,k,ispec)
             mu_buf(ipoint) = mustore(i,j,k,ispec)
             rho_buf(ipoint) = rhostore(i,j,k,ispec)

          enddo
       endif
    enddo

    close(IOUT)
    nbrec = sum(nb_gll_per_proc)
    ! dimensions for hdf5
    coords_dimsf(1) = 3
    coords_dimsf(2) = nbrec
    params_dimsf(1) = nbrec
    coords_dimsm(1) = 3
    coords_dimsm(2) = nb_gll_per_proc(myrank+1)
    params_dimsm(1) = nb_gll_per_proc(myrank+1)
    nbrec_dims = (/1/)
    spec_dims = (/sizeprocs/)
    coords_countf = (/3, nb_gll_per_proc(myrank+1)/)
    params_countf = (/nb_gll_per_proc(myrank+1)/)
    coords_offsetf = (/0, offset_per_proc(myrank+1)/)
    params_offsetf = (/offset_per_proc(myrank+1)/)


    call synchronize_all()
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS) then
      !! MPC write a hdf5 file for Instaseis input
      ! Initialize hdf5 interface
      call h5open_f(error)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio(plist_id)

      ! Open the file collectively.
      call h5fopen_f(hdf5_file, H5F_ACC_RDWR_F, file_id, error, &
                    access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      ! Open existing group
      call h5gopen_f(file_id, grp_local, grp_local_id, error)
      ! Create new group
      call h5gcreate_f(file_id, grp_params, grp_params_id, error)

      ! Create data spaces for the datasets, the dims are the entire arrays
      call h5screate_simple_f(coords_rank, coords_dimsf, dspace_coords_id, error)
      call h5screate_simple_f(coords_rank, coords_dimsf, dspace_normals_id, error)
      call h5screate_simple_f(params_rank, params_dimsf, dspace_kappa_id, error)
      call h5screate_simple_f(params_rank, params_dimsf, dspace_mu_id, error)
      call h5screate_simple_f(params_rank, params_dimsf, dspace_rho_id, error)
      call h5screate_simple_f(nbrec_rank, nbrec_dims, aspace_nbrec_id, error)
      call h5screate_simple_f(spec_rank, spec_dims, aspace_spec_id, error)
      call h5screate_simple_f(spec_rank, spec_dims, aspace_spec2_id, error)

      ! Create datasets with default properties.
      call h5dcreate_f(grp_local_id, dset_coords, H5T_NATIVE_REAL, &
                       dspace_coords_id, dset_coords_id, error)
      call h5dcreate_f(grp_local_id, dset_normals, H5T_NATIVE_REAL, &
                       dspace_normals_id, dset_normals_id, error)
      call h5dcreate_f(grp_params_id, dset_kappa, H5T_NATIVE_REAL, &
                       dspace_kappa_id, dset_kappa_id, error)
      call h5dcreate_f(grp_params_id, dset_mu, H5T_NATIVE_REAL, &
                       dspace_mu_id, dset_mu_id, error)
      call h5dcreate_f(grp_params_id, dset_rho, H5T_NATIVE_REAL, &
                       dspace_rho_id, dset_rho_id, error)
      call h5acreate_f(grp_local_id, attr_nbrec, H5T_NATIVE_INTEGER, &
                       aspace_nbrec_id, attr_nbrec_id, error)
      call h5acreate_f(grp_local_id, attr_spec, H5T_NATIVE_INTEGER, &
                      aspace_spec_id, attr_spec_id, error)
      call h5acreate_f(grp_local_id, attr_spec2, H5T_NATIVE_INTEGER, &
                      aspace_spec2_id, attr_spec2_id, error)

      call h5sselect_hyperslab_f(dspace_coords_id, H5S_SELECT_SET_F, &
                                 coords_offsetf, coords_countf, error)
      call h5sselect_hyperslab_f(dspace_normals_id, H5S_SELECT_SET_F, &
                                 coords_offsetf, coords_countf, error)
      call h5sselect_hyperslab_f(dspace_kappa_id, H5S_SELECT_SET_F, &
                                 params_offsetf, params_countf, error)
      call h5sselect_hyperslab_f(dspace_mu_id, H5S_SELECT_SET_F, &
                                 params_offsetf, params_countf, error)
      call h5sselect_hyperslab_f(dspace_rho_id, H5S_SELECT_SET_F, &
                                 params_offsetf, params_countf, error)

      ! Create memory dataspace, the dims are the arrays per proc
      call h5screate_simple_f(coords_rank, coords_dimsm, mspace_coords_id, error)
      call h5screate_simple_f(coords_rank, coords_dimsm, mspace_normals_id, error)
      call h5screate_simple_f(params_rank, params_dimsm, mspace_kappa_id, error)
      call h5screate_simple_f(params_rank, params_dimsm, mspace_mu_id, error)
      call h5screate_simple_f(params_rank, params_dimsm, mspace_rho_id, error)

      ! Create property list for independent dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio(plist_id)

      ! Write data and attributes
      call h5dwrite_f(dset_coords_id, H5T_NATIVE_REAL, coords_buf, &
          coords_dimsm, error, mspace_coords_id, dspace_coords_id, &
          xfer_prp = plist_id)
      call h5dwrite_f(dset_normals_id, H5T_NATIVE_REAL, normals_buf, &
          coords_dimsm, error, mspace_normals_id, dspace_normals_id, &
          xfer_prp = plist_id)
      call h5dwrite_f(dset_kappa_id, H5T_NATIVE_REAL, kappa_buf, params_dimsm, &
          error, mspace_kappa_id, dspace_kappa_id, xfer_prp = plist_id)
      call h5dwrite_f(dset_mu_id, H5T_NATIVE_REAL, mu_buf, params_dimsm, &
          error, mspace_mu_id, dspace_mu_id, xfer_prp = plist_id)
      call h5dwrite_f(dset_rho_id, H5T_NATIVE_REAL, rho_buf, params_dimsm, &
          error, mspace_rho_id, dspace_rho_id, xfer_prp = plist_id)
      call h5awrite_f(attr_nbrec_id, H5T_NATIVE_INTEGER, nbrec, nbrec_dims, error)
      call h5awrite_f(attr_spec_id, H5T_NATIVE_INTEGER, nb_gll_per_proc, &
                      spec_dims, error)
      call h5awrite_f(attr_spec2_id, H5T_NATIVE_INTEGER, offset_per_proc, &
                      spec_dims, error)

      ! Close hdf5 groups, datasets, dataspaces attributes...
      call h5sclose_f(dspace_coords_id, error)
      call h5sclose_f(dspace_normals_id, error)
      call h5sclose_f(dspace_kappa_id, error)
      call h5sclose_f(dspace_mu_id, error)
      call h5sclose_f(dspace_rho_id, error)
      call h5sclose_f(mspace_coords_id, error)
      call h5sclose_f(mspace_normals_id, error)
      call h5sclose_f(mspace_kappa_id, error)
      call h5sclose_f(mspace_mu_id, error)
      call h5sclose_f(mspace_rho_id, error)
      call h5sclose_f(aspace_nbrec_id, error)
      call h5sclose_f(aspace_spec_id, error)
      call h5sclose_f(aspace_spec2_id, error)
      call h5dclose_f(dset_coords_id, error)
      call h5dclose_f(dset_normals_id, error)
      call h5dclose_f(dset_kappa_id, error)
      call h5dclose_f(dset_mu_id, error)
      call h5dclose_f(dset_rho_id, error)
      call h5aclose_f(attr_nbrec_id, error)
      call h5aclose_f(attr_spec_id, error)
      call h5aclose_f(attr_spec2_id, error)
      call h5gclose_f(grp_local_id, error)
      call h5gclose_f(grp_params_id, error)
      call h5pclose_f(plist_id, error)
      ! Close the file.
      call h5fclose_f(file_id, error)

      ntime = NSTEP

      !! MPC Create a file that is later going to be filled out by specfem
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio(plist_id)
      ! Create the file collectively.
      call h5fcreate_f(spec_file, H5F_ACC_TRUNC_F, spec_file_id, error, &
                       access_prp = plist_id)

      call h5pclose_f(plist_id, error)
      ! Create new group
      call h5gcreate_f(spec_file_id, grp_coords, spec_grp_id, error)

      d_dimsf(1) = 3
      d_dimsf(2) = ntime
      d_dimsf(3) = nbrec
      s_dimsf(1) = 6
      s_dimsf(2) = ntime
      s_dimsf(3) = nbrec

      ! Create data spaces for the datasets, the dims are the entire arrays
      call h5screate_simple_f(v_rank, d_dimsf, dspace_spec_d_id, error)
      call h5screate_simple_f(s_rank, s_dimsf, dspace_spec_s_id, error)
      ! Create datasets with default properties.
      call h5dcreate_f(spec_grp_id, dset_spec_d, H5T_NATIVE_REAL, &
                       dspace_spec_d_id, dset_spec_d_id, error)
      call h5dcreate_f(spec_grp_id, dset_spec_s, H5T_NATIVE_REAL, &
                       dspace_spec_s_id, dset_spec_s_id, error)

      call h5dclose_f(dset_spec_d_id, error)
      call h5dclose_f(dset_spec_s_id, error)
      call h5gclose_f(spec_grp_id, error)
      ! Close the file.
      call h5fclose_f(spec_file_id, error)

      ! Close hdf5 interface
      call h5close_f(error)
    endif ! (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_INSTASEIS)
  endif ! (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)

  end subroutine save_arrays_solver_files
