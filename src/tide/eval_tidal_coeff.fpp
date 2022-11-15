! Program  : eval_tidal_coeff
! Purpose  : evaluate tidal coupling coefficients
!
! Copyright 2020-2022 Rich Townsend & The GYRE Team
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

  use gyre_hom_model
  use gyre_math
  use gyre_model_par
  use gyre_orbit_par
  use gyre_func
  use gyre_tidal_coeff

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

  type(model_par_t) :: ml_p
  type(hom_model_t) :: ml
  type(orbit_par_t) :: or_p
  real(WP)          :: q
  real(WP)          :: Omega_orb
  real(WP)          :: f

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

  ! Set up model and orbit parameters

  ml_p = model_par_t()

  ml = hom_model_t(ml_p)

  q = 1._WP
  Omega_orb = sqrt(R_a**3*(1._WP + q))

  or_p = orbit_par_t(Omega_orb=Omega_orb, q=q, e=e)

  ! Evaluate the coefficient

  select case (coeff)
  case ('cbar')
     f = tidal_cbar(ml, or_p, l, m, k)
  case ('Gbar_1')
     f = tidal_Gbar_1(ml, or_p, l, m, k)
  case ('Gbar_2')
     f = tidal_Gbar_2(ml, or_p, l, m, k)
  case ('Gbar_3')
     f = tidal_Gbar_3(ml, or_p, l, m, k)
  case ('Gbar_4')
     f = tidal_Gbar_4(ml, or_p, l, m, k)
  case ('X')
     f = hansen_X(or_p, l, m, k)
  case ('X_QP')
     f = hansen_X_QP(or_p, l, m, k)
  case ('Y')
     f = REAL(spherical_Y(l, m, HALFPI, 0._WP))
  case ('Y*')
     f = REAL(CONJG(spherical_Y(l, m, HALFPI, 0._WP)))
  case default
     $ABORT(Invalid coeff; should be one of [c|G_1|G_2|G_3|G_4|X|Y|Y*])
  end select

  print *, f

end program eval_tidal_coeff
