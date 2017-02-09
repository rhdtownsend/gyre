! Program  : poly_to_fgong
! Purpose  : convert a polytrope to FGONG format
!
! Copyright 2015-2017 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

$include 'core.inc'

program poly_to_fgong

  ! Uses

  use core_kinds
  use core_system

  use gyre_constants
  use gyre_grid
  use gyre_model
  use gyre_model_par
  use gyre_poly_file
  use gyre_poly_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: IVERS = 1300
  integer, parameter :: ICONST = 15
  integer, parameter :: IVAR = 40

  ! Variables

  character(:), allocatable :: in_filename
  character(:), allocatable :: out_filename
  logical                   :: drop_outer

  type(model_par_t)       :: ml_p
  class(model_t), pointer :: ml
  type(grid_t)            :: gr
  real(WP), allocatable   :: V_2(:)
  real(WP), allocatable   :: As(:)
  real(WP), allocatable   :: U(:)
  real(WP), allocatable   :: c_1(:)
  real(WP), allocatable   :: Gamma_1(:)
  real(WP), allocatable   :: M_r(:)
  real(WP), allocatable   :: P(:)
  real(WP), allocatable   :: rho(:)
  integer                 :: n_k
  real(WP), allocatable   :: glob(:)
  real(WP), allocatable   :: var(:,:)
  integer                 :: unit
  integer                 :: k

  ! Read parameters

  $ASSERT(n_arg() == 3,Syntax: poly_to_fgong in_filename out_filename drop_outer)

  call get_arg(1, in_filename)
  call get_arg(2, out_filename)
  call get_arg(3, drop_outer)

  ! Read the polytrope data

  ml_p%file = TRIM(in_filename)

  call read_poly_model(ml_p, ml)

  ! Set up the grid

  gr = ml%grid()

  if (drop_outer) then
     gr = grid_t(gr%pt(:gr%n_k-1)%x)
  endif

  ! Extract data from the model

  associate (pt => gr%pt)

    ! Dimensionless structure variables

    do k = 1, gr%n_k
       V_2 = ml%coeff(I_V_2, pt(k))
       As = ml%coeff(I_AS, pt(k))
       U = ml%coeff(I_U, pt(k))
       c_1 = ml%coeff(I_C_1, pt(k))
       Gamma_1 = ml%coeff(I_GAMMA_1, pt(k))
    end do

    ! Physical structure variables

    $if ($GFORTRAN_PR_49636)
    M_r = M_SUN*(gr%pt%x**3/c_1)
    $else
    M_r = M_SUN*(pt%x**3/c_1)
    $endif

    P = (G_GRAVITY*M_SUN**2/(4._WP*PI*R_SUN**4))*&
        (U/(c_1**2*V_2))

    rho = (M_SUN/(4._WP*PI*R_SUN**3))*(U/c_1)

  end associate

  ! Store into var array

  n_k = gr%n_k

  allocate(glob(ICONST))
  allocate(var(IVAR,n_k))

  glob(1) = M_SUN
  glob(2) = R_SUN
  glob(3:) = 0._WP

  var(1,:) = gr%pt%x*R_SUN

  where (M_r /= 0._WP)
     var(2,:) = LOG(M_r/M_SUN)
  elsewhere
     var(2,:) = LOG(1E-38_WP)
  endwhere

  var(3,:) = 0._WP
  var(4,:) = P
  var(5,:) = rho
  var(6:9,:) = 0._WP
  var(10,:) = Gamma_1
  var(11:14,:) = 0._WP
  var(15,:) = As
  var(16,:) = 0._WP

  var = var(:,n_k:1:-1)

  ! Write out the FGONG file

  open(NEWUNIT=unit, FILE=out_filename, STATUS='REPLACE')

  write(unit, 100) 'Fee'
  write(unit, 100) 'Fi'
  write(unit, 100) 'Fo'
  write(unit, 100) 'Fum'
100 format(A)

  write(unit, 110) n_k, ICONST, IVAR, IVERS
110 format(4I10)

  write(unit, 120) glob
120  format(1P,5(1X,E26.18E3))

  do k = 1, n_k
     write(unit, 120) var(:,k)
  end do

  close(unit)

  ! Finish

end program poly_to_fgong
