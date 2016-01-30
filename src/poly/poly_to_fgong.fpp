! Program  : poly_to_fgong
! Purpose  : convert a polytrope to FGONG format
!
! Copyright 2015-2016 Rich Townsend
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

  type(model_par_t)         :: ml_p
  class(model_t), pointer   :: ml
  integer, allocatable      :: s(:)
  real(WP), allocatable     :: x(:)
  real(WP), allocatable     :: V_2(:)
  real(WP), allocatable     :: As(:)
  real(WP), allocatable     :: U(:)
  real(WP), allocatable     :: c_1(:)
  real(WP), allocatable     :: Gamma_1(:)
  real(WP), allocatable     :: M_r(:)
  real(WP), allocatable     :: P(:)
  real(WP), allocatable     :: rho(:)
  integer                   :: n
  real(WP), allocatable     :: glob(:)
  real(WP), allocatable     :: var(:,:)
  integer                   :: unit
  integer                   :: i

  ! Read parameters

  $ASSERT(n_arg() == 3,Syntax: poly_to_fgong in_filename out_filename drop_outer)

  call get_arg(1, in_filename)
  call get_arg(2, out_filename)
  call get_arg(3, drop_outer)

  ! Read the polytrope data

  ml_p%file = TRIM(in_filename)

  call read_poly_model(ml_p, ml)

  ! Set up the grid

  call set_grid(ml, drop_outer, s, x)

  ! Calculate dimensionless structure data

  V_2 = ml%V_2(s, x)
  As = ml%As(s, x)
  U = ml%U(s, x)
  c_1 = ml%c_1(s, x)
  Gamma_1 = ml%Gamma_1(s, x)

  ! Calculate physical data

  M_r = M_SUN*(x**3/c_1)

  P = (G_GRAVITY*M_SUN**2/(4._WP*PI*R_SUN**4))*&
       (U/(c_1**2*V_2))

  rho = (M_SUN/(4._WP*PI*R_SUN**3))*(U/c_1)

  ! Store into var array

  n = SIZE(s)

  allocate(glob(ICONST))
  allocate(var(IVAR,n))

  glob(1) = M_SUN
  glob(2) = R_SUN
  glob(3:) = 0._WP

  var(1,:) = x*R_SUN

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

  var = var(:,n:1:-1)

  ! Write out the FGONG file

  open(NEWUNIT=unit, FILE=out_filename, STATUS='REPLACE')

  write(unit, 100) 'Fee'
  write(unit, 100) 'Fi'
  write(unit, 100) 'Fo'
  write(unit, 100) 'Fum'
100 format(A)

  write(unit, 110) n, ICONST, IVAR, IVERS
110 format(4I10)

  write(unit, 120) glob
120  format(1P,5(X,E26.18E3))

  do i = 1, n
     write(unit, 120) var(:,i)
  end do

  close(unit)

  ! Finish

contains

  subroutine set_grid (ml, drop_outer, s, x)

    class(model_t), pointer, intent(in) :: ml
    logical, intent(in)                 :: drop_outer
    integer, allocatable, intent(out)   :: s(:)
    real(WP), allocatable, intent(out)  :: x(:)

    integer               :: s_
    real(WP), allocatable :: x_s(:)

    ! Set up the grid

    allocate(s(0))
    allocate(x(0))

    seg_loop : do s_ = 1, ml%n_s

       x_s = ml%x_base(s_)

       s = [s,SPREAD(s_, 1, SIZE(x_s))]
       x = [x,x_s]

    end do seg_loop

    if (drop_outer) then
       s = s(:SIZE(s)-1)
       x = x(:SIZE(x)-1)
    endif

    ! Finish

    return

  end subroutine set_grid
    
end program poly_to_fgong
