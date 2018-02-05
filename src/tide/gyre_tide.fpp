! Program  : gyre_tide
! Purpose  : tide-related functions
!
! Copyright 2018 Rich Townsend
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

module gyre_tide

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: c_lmk

  ! Procedures

contains

  function c_lmk (r_a, ec, l, m, k)

    real(WP), intent(in) :: r_a
    real(WP), intent(in) :: ec
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: c_lmk

    integer, parameter :: N = 1024

    integer  :: i
    real(WP) :: ua(N)
    real(WP) :: Ea
    real(WP) :: Ma
    real(WP) :: y(N)
    real(WP) :: X
    real(WP) :: P
    real(WP) :: F
    integer  :: j

    ! Evaluate the tidal Fourier coefficient c_lmk (e.g., Smeyers,
    ! Willems & Van Hoolst 1998, A&A, 335, 622)

    if (ABS(m) > l) then

       c_lmk = 0._WP

    else

       if (MOD(l-ABS(m), 2) == 0) then

          ! Calculate the Hansen coefficient via eqn. (22)

          do i = 1, N

             ua(i) = (i-1)*PI/(N-1)

             Ea = 2._WP*ATAN(SQRT((1._WP-ec)/(1._WP+ec))*TAN(ua(i)/2._WP))

             Ma = Ea - ec*SIN(Ea)

             y(i) = (1._WP + ec*COS(ua(i)))**(l-1) * COS(k*Ma + m*ua(i))

          end do

          X = integrate(ua, y)/(1._WP - ec**2)**(l-0.5_WP)/PI

          ! Calculate the associated Legendre function P|m|l(0) and
          ! the factorial term F = (l-|m|)!/(l+|m|)!

          P = 1._WP
          F = 1._WP

          do j = 1, l+ABS(m)

             if (MOD(j, 2) == 0) then
                if (j <= l-ABS(m)) P = P/j
             else
                if (j <= l+ABS(m)-1) P = P*j
             endif

             if (j <= l+ABS(m)) F = F/j
             if (j <= l-ABS(m)) F = F*j

          end do

          P = (-1)**((l+ABS(m))/2)*P

          ! Calculate the Fourier coefficient via eqn. (18)

          c_lmk = F*P*R_a**(l-2)*X

       else

          c_lmk = 0._WP

       endif

    end if

    ! Finish

    return

  end function c_lmk

end module gyre_tide
