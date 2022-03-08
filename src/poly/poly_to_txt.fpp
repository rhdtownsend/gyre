! Program  : poly_to_txt
! Purpose  : convert a polytrope to a simple txt format
!
! Copyright 2019 Rich Townsend & The GYRE Team
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

program poly_to_txt

  ! Uses

  use core_kinds
  use core_system

  use gyre_grid
  use gyre_model
  use gyre_model_par
  use gyre_poly_file
  use gyre_poly_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

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
  integer                 :: n_p
  integer                 :: unit
  integer                 :: p

  ! Read parameters

  $ASSERT(n_arg() == 3,Syntax: poly_to_txt in_filename out_filename drop_outer)

  call get_arg(1, in_filename)
  call get_arg(2, out_filename)
  call get_arg(3, drop_outer)

  ! Read the polytrope data

  ml_p%file = TRIM(in_filename)

  call read_poly_model(ml_p, ml)

  ! Set up the grid

  gr = ml%grid()

  if (drop_outer) then
     gr = grid_t(gr%pt(:gr%n_p-1)%x)
  endif

  ! Extract data from the model

  ! Dimensionless structure variables

  allocate(V_2(gr%n_p))
  allocate(As(gr%n_p))
  allocate(U(gr%n_p))
  allocate(c_1(gr%n_p))
  allocate(Gamma_1(gr%n_p))

  do p = 1, gr%n_p
     associate (pt => gr%pt(p))
       V_2(p) = ml%coeff(I_V_2, pt)
       As(p) = ml%coeff(I_AS, pt)
       U(p) = ml%coeff(I_U, pt)
       c_1(p) = ml%coeff(I_C_1, pt)
       Gamma_1(p) = ml%coeff(I_GAMMA_1, pt)
     end associate
  end do

  ! Write out the txt file

  open(NEWUNIT=unit, FILE=out_filename, STATUS='REPLACE')

  write(unit, 100) 'x V_2 As U c_1 Gamma_1'
100 format(A)

  do p = 1, gr%n_p
     write(unit, 110) gr%pt(p)%x, V_2(p), As(p), U(p), c_1(p), Gamma_1(p)
110  format(6(1X,E26.18E3))
  end do

  close(unit)

  ! Finish

end program poly_to_txt
