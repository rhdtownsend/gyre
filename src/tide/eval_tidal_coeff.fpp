! Program  : eval_tidal_coeff
! Purpose  : evaluate tidal coupling coefficients
!
! Copyright 2020 Rich Townsend & The GYRE Team
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

program eval_tidal_coeff

  ! Uses

  use core_kinds
  use core_system

  use gyre_math
  use gyre_func
  use gyre_tide_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable :: coeff
  real(WP)                  :: R_a
  real(WP)                  :: e
  integer                   :: l
  integer                   :: m
  integer                   :: k

  real(WP) :: f

  ! Get arguments

  if (n_arg() /= 6) stop '** syntax: eval_tidal_coeff coeff R_a e [n|l] m k'

  call get_arg(1, coeff)
  call get_arg(2, R_a)
  call get_arg(3, e)
  call get_arg(4, l)
  call get_arg(5, m)
  call get_arg(6, k)

  ! Initialize

  call init_math()

  ! Evaluate the coefficient

  select case (coeff)
  case ('c')
     f = tidal_c(R_a, e, l, m, k)
  case ('G_1')
     f = secular_G_1(R_a, e, l, m, k)
  case ('G_2')
     f = secular_G_2(R_a, e, l, m, k)
  case ('G_3')
     f = secular_G_3(R_a, e, l, m, k)
  case ('G_4')
     f = secular_G_4(R_a, e, l, m, k)
  case ('X')
     f = hansen_X(e, l, m, k)
  case ('Y')
     f = REAL(spherical_Y(l, m, HALFPI, 0._WP))
  case ('Y*')
     f = REAL(CONJG(spherical_Y(l, m, HALFPI, 0._WP)))
  case default
     $ABORT(Invalid coeff; should be one of [c|G_1|G_2|G_3|G_4|X|Y|Y*])
  end select

  print *, f

end program eval_tidal_coeff
  
