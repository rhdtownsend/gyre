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

  ! Parameters

  real(WP), parameter :: TOL = 1E-12_WP

  ! Access specifiers

  private

  public :: tidal_c
  public :: tidal_kappa
  public :: hansen_X
  public :: hansen_X_tilde
  public :: hansen_X_hat

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

    c = (4._WP*PI/(2*l+1))*R_a**(l-2)*REAL(CONJG(spherical_Y(l, m, HALFPI, 0._WP)))*hansen_X(e, -(l+1), -m, k)

    ! Finish

    return

  end function tidal_c

  !****

  function tidal_kappa (l, m, k) result (kappa)

    integer, intent(in) :: l
    integer, intent(in) :: m
    integer, intent(in) :: k
    real(WP)            :: kappa

    $ASSERT(k >= 0,Invalid k)

    ! Evaluate the kappa_lmk function, using the expression given in
    ! equation (53) of Willems et al. (2010)

    if (k == 0) then

       if (m == 0) then
          kappa = 0.5_WP
       elseif (m >= 1) then
          kappa = 1._WP
       else
          kappa = 0._WP
       endif

    else

       kappa = 1._WP

    endif

    ! Finish

    return

  end function tidal_kappa

  !****

  function hansen_X (e, n, m, k) result (X)

    real(WP), intent(in) :: e
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: X

    real(WP) :: A

    $ASSERT(e >= 0,Invalid e)
    $ASSERT(e < 1,Invalid e)

    ! Evaluate the Hansen coefficient X_nmk

    if (e**2 > EPSILON(0._WP)) then

       A = (1._WP - e**2)**(n+1.5_WP)/TWOPI

       X = orbit_quad(hansen_X_f_, e, TOL)

    else

       if (m == k) then
          X = 1._WP
       else
          X = 0._WP
       endif

    endif

    ! Finish

    return

  contains

    function hansen_X_f_ (ua, Ma, e) result (f)

      real(WP), intent(in) :: ua
      real(WP), intent(in) :: Ma
      real(WP), intent(in) :: e
      real(WP)             :: f

      ! Set up the integrand

      f = A*COS(m*ua - k*Ma)/(1._WP + e*COS(ua))**(n+2)

      ! Finish

      return

    end function hansen_X_f_

  end function hansen_X

  !****

  function hansen_X_tilde (e, n, m, k) result (X_tilde)

    real(WP), intent(in) :: e
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: X_tilde

    real(WP) :: A

    $ASSERT(e >= 0,Invalid e)
    $ASSERT(e < 1,Invalid e)

    ! Evaluate the tilded Hansen coefficient \tilde{X}_nmk

    if (e**2 > EPSILON(0._WP)) then

       A = (1._WP - e**2)**(n+1.5_WP)/TWOPI

       X_tilde = orbit_quad(hansen_X_tilde_f_, e, TOL)

    else

       if (ABS(m - k) == 1) then
          X_tilde = 0.5_WP*SIGN(1, m-k)
       else
          X_tilde = 0._WP
       endif

    endif

    ! Finish

    return

  contains

    function hansen_X_tilde_f_ (ua, Ma, e) result (f)

      real(WP), intent(in) :: ua
      real(WP), intent(in) :: Ma
      real(WP), intent(in) :: e
      real(WP)             :: f

      ! Set up the integrand

      f = A*SIN(m*ua - k*Ma)*SIN(ua)/(1._WP + e*COS(ua))**(n+2)

      ! Finish

      return

    end function hansen_X_tilde_f_

  end function hansen_X_tilde

  !****

  function hansen_X_hat (e, n, m, k) result (X_hat)

    real(WP), intent(in) :: e
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: X_hat

    real(WP) :: A

    $ASSERT(e >= 0,Invalid e)
    $ASSERT(e < 1,Invalid e)

    ! Evaluate the hatted Hansen coefficient \hat{X}_nmk

    if (e**2 > EPSILON(0._WP)) then

       A = (1._WP - e**2)**(n+1.5_WP)/TWOPI

       X_hat = orbit_quad(hansen_X_hat_f_, e, TOL)

    else

       if (ABS(m - k) == 1) then
          X_hat = 0.5_WP
       else
          X_hat = 0._WP
       endif

    endif

    ! Finish

    return

  contains

    function hansen_X_hat_f_ (ua, Ma, e) result (f)

      real(WP), intent(in) :: ua
      real(WP), intent(in) :: Ma
      real(WP), intent(in) :: e
      real(WP)             :: f

      ! Set up the integrand

      f = A*COS(m*ua - k*Ma)*COS(ua)/(1._WP + e*COS(ua))**(n+2)

      ! Finish

      return

    end function hansen_X_hat_f_

  end function hansen_X_hat

  !****

  function orbit_quad (f, e, tol) result (I)

    interface
       function f (ua, Ma, e)
         use core_kinds
         real(WP), intent(in) :: ua
         real(WP), intent(in) :: Ma
         real(WP), intent(in) :: e
         real(WP)             :: f
       end function f
    end interface
    real(WP), intent(in) :: e
    real(WP), intent(in) :: tol
    real(WP)             :: I

    integer, parameter :: N_MAX = 2**20

    integer  :: N
    real(WP) :: e_fac
    real(WP) :: S
    integer  :: j
    real(WP) :: ua
    real(WP) :: Ea
    real(WP) :: Ma
    real(WP) :: dI
    real(WP) :: dI_max

    ! Perform "orbital quadrature" on the function f(ua, Ma, e, l, m,
    ! k). This involves calculating the integral I = int(f(ua),
    ! ua->0,2pi), where ua is the true anomaly.  Along with ua, the
    ! corresponding mean anomaly Ma(ua) is passed into f

    ! Set up initial values

    N = 2

    I = PI*(f(0._WP, 0._WP, e) + f(PI, PI, e))

    dI_max = 0._WP

    ! Refine the quadrature by repeatedly doubling N

    e_fac = SQRT((1._WP - e)*(1._WP + e))

    refine_loop : do

       ! Update the value of the quadrature

       S = 0._WP

       update_loop : do j = 1, N

          ! Evaluate the true anomaly

          ua = TWOPI*(j-0.5_WP)/N

          ! Evaluate the eccentric anomaly

          Ea = ATAN2(e_fac*SIN(ua), e + COS(ua))

          ! Evaluate the mean anomaly

          Ma = Ea - e*SIN(Ea)

          ! Add the contribution

          S = S + f(ua, Ma, e)

       end do update_loop

       dI = PI*S/N - 0.5_WP*I
       I = I + dI

       N = 2*N

       ! Check for convergence; the current correction is smaller than
       ! tol times the largest

       $ASSERT(N <= N_MAX,Too many iterations)

       dI_max = MAX(dI_max, ABS(dI))

       if (ABS(dI) < tol*dI_max) exit refine_loop

    end do refine_loop

    ! Finish

    return

  end function orbit_quad

end module gyre_tide_util
