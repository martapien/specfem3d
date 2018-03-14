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

!--------------------------------------------------------------------------------------------------
!
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------
!
!
! example model to couple with an injection technique (FK,DSM,AXISEM)
!

  module model_coupled_par

  ! VM VM my model for DSM coupling
  use constants

  ! VM VM
  double precision, dimension (:,:), allocatable :: vpv_1D,vsv_1D,density_1D
  double precision, dimension (:), allocatable :: zlayer
  double precision, dimension (:), allocatable :: smooth_vp,smooth_vs
  integer :: ilayer,nlayer,ncoeff,ndeg_poly
  double precision :: ZREF,OLON,OLAT

  !! gaussain perturbation definition
  double precision :: sx, sy, sz
  double precision :: Ampl_pert_vp_read, Ampl_pert_vs_read, Ampl_pert_rho_read
  double precision :: Ampl_pert_vp, Ampl_pert_vs, Ampl_pert_rho
  double precision :: x_center_gauss, y_center_gauss, z_center_gauss


  end module model_coupled_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_coupled_broadcast(myrank)

! standard routine to setup model

  use model_coupled_par

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH, ANISOTROPY, ATTENUATION

  implicit none

  integer :: myrank

  ! local parameters
  integer :: idummy

  ! dummy to ignore compiler warnings
  idummy = myrank

  ! safety check
  if (.not. (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)) then
    print *,'Error: model coupling requires coupling with injection technique or mesh a chunk of the earth'
    stop 'Error model coupling'
  endif

  call read_model_for_coupling_or_chunk()

  end subroutine model_coupled_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_for_coupling_or_chunk()

  use model_coupled_par !! VM VM custom subroutine for coupling with DSM
  use shared_parameters, only:  ANISOTROPY, ATTENUATION, ADD_GAUSSIAN_PERT_ABSOLUTE

  implicit none

  character(len=256) :: filename
  character(len=10)  :: line
  character(len=250) :: model1D_file
  integer            :: i
  double precision   :: ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD
  double precision   :: lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision   :: radius_of_box_top
  logical            :: buried_box
  integer            :: nel_lat, nel_lon, nel_depth

  !! reading chunk parameters
  open(27,file='MESH/ParFileMeshChunk',action='read')
  read(27,'(a)') line
  read(27,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  read(27,'(a)') line
  read(27,*) lon_center_chunk, lat_center_chunk, chunk_azi
  read(27,'(a)') line
  read(27,*) chunk_depth
  read(27,'(a)') line
  read(27,*) nel_lon,nel_lat, nel_depth
  read(27,'(a)') line
  read(27,'(a)') model1D_file
  read(27,'(a)') line
  read(27,*) buried_box
  if (buried_box) then
     read(27,'(a)') line
     read(27,*) radius_of_box_top
     radius_of_box_top =  radius_of_box_top * 1000.
  else
     radius_of_box_top = 6371000.
  endif
  close(27)

  OLON = lon_center_chunk
  OLAT = lat_center_chunk
  ZREF = radius_of_box_top

  !write(*,*) " Reading 1D model "

  !! reading 1D reference model given by polynomial coeffiscients
  if (ANISOTROPY) then

     !! read aniso model
     write(*,*) " external 1D anisotropic model not defined yet "
     stop

     if (ATTENUATION) then
        !! TO DO  ...
        !! read 1D model with attenuation
        write(*,*) " external 1D vsico-elastic model not defined yet "
        stop

     else
        !! TO DO  ...
     end if

  else

     if (ATTENUATION) then

        !! TO DO  ...
        !! read 1D model with attenuation
        write(*,*) " external 1D vsico-elastic model not defined yet "
        stop

     else

        !! isotropic 1D model:
        filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//trim(model1D_file)
        open(27,file=trim(filename))
        read(27,*) nlayer,ncoeff
        allocate(vpv_1D(nlayer,ncoeff))
        allocate(vsv_1D(nlayer,ncoeff))
        allocate(density_1D(nlayer,ncoeff))
        allocate(zlayer(nlayer))
        do i=1,nlayer
           read(27,*) zlayer(i)
           read(27,*) vpv_1D(i,:)
           read(27,*) vsv_1D(i,:)
           read(27,*) density_1D(i,:)
        enddo
        close(27)

     end if

  end if


  if (ADD_GAUSSIAN_PERT_ABSOLUTE) then
    !! MPC reading gaussian pert parameters
    open(27,file='MESH/ParFileGaussianPert',action='read')
    read(27, '(a)') line
    read(27, '(a)') line
    read(27, *)     x_center_gauss, y_center_gauss, z_center_gauss
    read(27, '(a)') line
    read(27, *)     sx, sy, sz
    read(27, '(a)') line
    read(27, *)     Ampl_pert_vp_read
    read(27, '(a)') line
    read(27, *)     Ampl_pert_vs_read
    read(27, '(a)') line
    read(27, *)     Ampl_pert_rho_read
    close(27)

    Ampl_pert_vp  = dble(Ampl_pert_vp_read)
    Ampl_pert_vs  = dble(Ampl_pert_vs_read)
    Ampl_pert_rho = dble(Ampl_pert_rho_read)
  endif

  !! hardcoded gaussian pert (todo need to read an input file)
  !! x_center_gauss=0.
  !! y_center_gauss=0.
  !! z_center_gauss=-100000.
  !! sx=100000.
  !! sy=100000.
  !! sz=50000.
  !! Ampl_pert_vp=1000.d0 !! abs olute amplitude perturbation
  !! Ampl_pert_vs=400.d0
  ! Ampl_pert_rho=400.d0

  end subroutine read_model_for_coupling_or_chunk

!----------------------------------------------------------------
!! !! ================= VM VM custom subroutine for DSM coupling
!----------------------------------------------------------------

  subroutine model_coupled_FindLayer(x,y,z)

  use model_coupled_par

  implicit none
  integer il
  double precision radius
  double precision :: x,y,z

  radius =  dsqrt(x**2 + y**2 + (z+zref)**2) / 1000.d0
  il = 1
  do while (radius > zlayer(il) .and. il < nlayer)
     il = il + 1
  enddo
  il = il - 1
  ilayer = il

  end subroutine model_coupled_FindLayer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_coupled_values(xmesh,ymesh,zmesh,rho,vp,vs)

! given a GLL point, returns super-imposed velocity model values

  use constants, only: CUSTOM_REAL

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE, MESH_A_CHUNK_OF_THE_EARTH, ADD_GAUSSIAN_PERT_ABSOLUTE

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho

  ! local parameters
  double precision :: x,y,z
  double precision :: radius

  ! safety check
  if (.not. (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)) then
      print *,'Error: model coupling requires coupling with injection technique or mesh a chunk of the earth'
    stop 'Error model coupling'
  endif

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh

  call  model_1D_coupling(x,y,z,rho,vp,vs,radius)


  if (ADD_GAUSSIAN_PERT_ABSOLUTE) then
      call add_gaussian_pert(x, y, z, rho, vp, vs)
  endif

  end subroutine model_coupled_values

!----------------------------------------------------------------


  subroutine add_gaussian_pert(x, y, z, rho, vp, vs)

    use constants, only: CUSTOM_REAL
    use  model_coupled_par

    double precision,       intent(in)    :: x, y, z
    real(kind=CUSTOM_REAL), intent(inout) :: rho, vp, vs
    double precision                      :: gauss_value

    gauss_value = exp( -0.5d0 * (    ((x - x_center_gauss) / sx)**2 + &
                                     ((y - y_center_gauss) / sy)**2 + &
                                     ((z - z_center_gauss) / sz)**2) )

    vp = vp + Ampl_pert_vp * gauss_value
    vs = vs + Ampl_pert_vs * gauss_value
    rho = rho +  Ampl_pert_rho * gauss_value

  end subroutine add_gaussian_pert

!----------------------------------------------------------------

  subroutine model_1D_coupling(x_eval,y_eval,z_eval,rho_final,vp_final,vs_final,r1)

  use model_coupled_par
  implicit none

  double precision r1,radius
  double precision rho,vp,vs
  double precision x_eval,y_eval,z_eval
  real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final

  double precision, parameter :: Xtol = 1d-2

  double precision :: xmin_pert, xmax_pert, ymin_pert, ymax_pert, zmin_pert, zmax_pert


  !! perturbed zone
!!$  xmin_pert = -100000.
!!$  xmax_pert =  100000.
!!$  ymin_pert = -100000.
!!$  ymax_pert =  100000.
!!$  zmin_pert = -180000.
!!$  zmax_pert =  -80000.


  radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+zref)**2)

!!$  if ( x_eval < 10. .and. x_eval > -10. .and. y_eval < 10. .and. y_eval > -10.) then
!!$     write(*,*)  z_eval / 1000. , radius / 1000. , 6371 + z_eval/ 1000. , zref
!!$  end if

  radius = radius / 1000.d0
  r1=radius

  ! get vp,vs and rho
  radius = radius / zlayer(nlayer)
  vp = Interpol(vpv_1D,ilayer,radius,nlayer)
  vs = Interpol(vsv_1D,ilayer,radius,nlayer)
  rho = Interpol(density_1D,ilayer,radius,nlayer)

  vp_final = vp * 1000.d0
  vs_final = vs * 1000.d0
  rho_final = rho * 1000.d0

  !! add perturbation
!!$  if ( x_eval > xmin_pert .and.  x_eval < xmax_pert .and. &
!!$       y_eval > ymin_pert .and.  y_eval < ymax_pert .and. &
!!$       z_eval > zmin_pert .and.  z_eval < zmax_pert ) then
!!$
!!$     vp_final = 1.1*vp_final
!!$     vs_final = 1.1*vs_final
!!$     rho_final = 1.1*rho_final
!!$
!!$  end if


  contains

    function Interpol(v,i,x,nl)

    implicit none
    integer i,nl
    double precision Interpol,x,v(nl,4)

    Interpol = v(i,1)+x*(v(i,2)+x*(v(i,3)+x*v(i,4)))

    end function Interpol

  end subroutine model_1D_coupling
