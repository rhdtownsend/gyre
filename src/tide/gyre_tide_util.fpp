! Program  : gyre_tide_util
! Purpose  : tide-related utility functions
!
! Copyright 2018-2019 Rich Townsend & The GYRE Team
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

module gyre_tide_util

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_func
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: tidal_c
  public :: hansen_X

  ! Procedures

contains

  function tidal_c (R_a, e, l, m, k) result (c)

    real(WP), intent(in)         :: R_a
    real(WP), intent(in)         :: e
    integer, intent(in)          :: l
    integer, intent(in)          :: m
    integer, intent(in)          :: k
    real(WP)                     :: c

    ! Evaluate the tidal potential coefficient

    c = (4._WP*PI/(2*l+1))*R_a**(l-2)*CONJG(spherical_Y(l, m, HALFPI, 0._WP))*hansen_X(e, -(l+1), -m, k)

    ! Finish

    return

  end function tidal_c

  !****

  function hansen_X (e, l, m, k) result (X)

    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: X

    integer, parameter :: N = 1024

    integer  :: i
    real(WP) :: ua(N)
    real(WP) :: Ea
    real(WP) :: Ma
    real(WP) :: y(N)

    ! Evaluate the Hansen coefficient X_lmk, using the transformed
    ! integral expression given in eqn. (22) of Smeyers, Willems & Van
    ! Hoolst (1998, A&A, 335, 622) (this expression avoids having to
    ! solve Kepler's equation)

    ! Set up the integrand

    do i = 1, N

       ua(i) = (i-1)*PI/(N-1)

       Ea = 2._WP*ATAN(SQRT((1._WP-e)/(1._WP+e))*TAN(ua(i)/2._WP))

       Ma = Ea - e*SIN(Ea)

       y(i) = COS(k*Ma - m*ua(i))/(1._WP + e*COS(ua(i)))**(l+2)

    end do

    ! Do the integral

    X = integrate(ua, y)*(1._WP - e**2)**(l+1.5_WP)/PI

    ! Finish

    return

  end function hansen_X

  !****


  function kappa_lmk (l, m, k)

    integer, intent(in) :: l
    integer, intent(in) :: m
    integer, intent(in) :: k
    real(WP)            :: kappa_lmk

    ! Evaluate the kappa_lmk function, using the expression given in
    ! equation (53) of Willems et al. (2010)

    if (k == 0) then

       if (m == 0) then
          kappa_lmk = 0.5_WP
       elseif (m >= 1) then
          kappa_lmk = 1._WP
       else
          kappa_lmk = 0._WP
       endif

    else

       kappa_lmk = 1._WP

    endif

    ! Finish

    return

  end function kappa_lmk

end module gyre_tide_util
