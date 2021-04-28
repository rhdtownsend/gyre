! Program  : gyre_tide_util
! Purpose  : tide-related utility functions
!
! Copyright 2018-2020 Rich Townsend & The GYRE Team
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
  use gyre_math
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  real(WP), parameter :: TOL = 1E-12_WP

  ! Access specifiers

  private

  public :: R_a
  public :: tidal_c
  public :: secular_G_1
  public :: secular_G_2
  public :: secular_G_3
  public :: secular_G_4
  public :: hansen_X

  ! Procedures

contains

  function tidal_R_a (Omega_orb, q) result (R_a)

     real(WP), intent(in) :: R_a
     real(WP), intent(in) :: q

     ! Evaluate the ratio of the stellar radius to the semi-major axis

     R_a = (Omega_orb**2/(1._WP + q))**(1._WP/3._WP)

     ! Finish

     return

  end function tidal_R_a

  !****

  function tidal_c (R_a, e, l, m, k) result (c)

    real(WP), intent(in) :: R_a
    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: c

    ! Evaluate the tidal potential coefficient. This differs from the
    ! definition e.g. in eq. (5) of [Willems:2010aa], due to the
    ! different spherical harmonic definition

    c = (4._WP*PI/(2*l+1))*R_a**(l-2)*REAL(CONJG(spherical_Y(l, m, HALFPI, 0._WP)))*hansen_X(e, -(l+1), -m, -k)

    ! Finish

    return

  end function tidal_c

  !****

  function secular_G_1 (R_a, e, l, m, k) result (G_1)

    real(WP), intent(in) :: R_a
    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: G_1

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_1m1
    real(WP) :: X_1p1
    real(WP) :: X_2m1
    real(WP) :: X_2p1

    ! Evaluate the secular evolution coefficient

    c = tidal_c(R_a, e, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_1m1 = hansen_X(e, -(l+1), -m-1, -k)
    X_1p1 = hansen_X(e, -(l+1), -m+1, -k)

    X_2m1 = hansen_X(e, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(e, -(l+2), -m+1, -k)

    G_1 = c*Y* &
         (0.5_WP*(l+1)*(X_2m1 + X_2p1) + 0.5_WP*m*(X_2m1 - X_2p1) + &
          0.5_WP*m/(1._WP - e**2)*(X_1m1 - X_1p1))*sqrt(1._WP - e**2)/e

    ! Finish

    return

  end function secular_G_1
  
  !****

  function secular_G_2 (R_a, e, l, m, k) result (G_2)

    real(WP), intent(in) :: R_a
    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: G_2

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_2m1
    real(WP) :: X_2p1
    real(WP) :: X_3

    ! Evaluate the secular evolution coefficient

    c = tidal_c(R_a, e, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_2m1 = hansen_X(e, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(e, -(l+2), -m+1, -k)

    X_3 = hansen_X(e, -(l+3), -m, -k)

    G_2 = -2._WP*c*Y* &
         (0.5_WP*(l+1)*e*(X_2m1 - X_2p1) + m*(1._WP - e**2)*X_3)/sqrt(1._WP - e**2)

    ! Finish

    return

  end function secular_G_2

  !****

  function secular_G_3 (R_a, e, l, m, k) result (G_3)

    real(WP), intent(in) :: R_a
    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: G_3

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_1
    real(WP) :: X_2m1
    real(WP) :: X_2p1
    real(WP) :: X_3

    ! Evaluate the secular evolution coefficient

    c = tidal_c(R_a, e, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_1 = hansen_X(e, -(l+1), -m, -k)

    X_2m1 = hansen_X(e, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(e, -(l+2), -m+1, -k)

    X_3 = hansen_X(e, -(l+3), -m, -k)

    G_3 = -c*Y* &
         (0.5_WP*(l+1)*e*(X_2m1 - X_2p1) + m*(1._WP - e**2)*X_3 - m*X_1)*sqrt(1._WP - e**2)/e

    ! Finish

    return

  end function secular_G_3

  !****

  function secular_G_4 (R_a, e, l, m, k) result (G_4)

    real(WP), intent(in) :: R_a
    real(WP), intent(in) :: e
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: G_4

    real(WP) :: c

    ! Evaluate the secular evolution coefficient

    c = tidal_c(R_a, e, l, m, k)

    G_4 = m*(2*l+1)/(4*PI)*(R_a)**(-l+2)*c**2

    ! Finish

    return

  end function secular_G_4

  !****

  function hansen_X (e, n, m, k) result (X)

    real(WP), intent(in) :: e
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: X

    real(WP) :: A
    real(WP) :: S

    $ASSERT(e >= 0,Invalid e)
    $ASSERT(e < 1,Invalid e)

    ! Evaluate the Hansen coefficient X_nmk

    if (e**2 > EPSILON(0._WP)) then

       A = (1._WP - e**2)**(n+1.5_WP)/TWOPI
       S = MAX((1._WP - e)**(-n-2), (1._WP + e)**(-n-2))

       X = A*S*orbit_quad(hansen_X_f_, e, TOL)

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

      f = cos(m*ua - k*Ma)/(S*(1._WP + e*cos(ua))**(n+2))

      ! Finish

      return

    end function hansen_X_f_

  end function hansen_X

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

    integer, parameter :: N_MAX = 2**24

    integer  :: N
    real(WP) :: e_fac
    real(WP) :: S
    integer  :: j
    real(WP) :: ua
    real(WP) :: Ea
    real(WP) :: Ma
    real(WP) :: dI
    real(WP) :: dI_conv

    ! Perform "orbital quadrature" on the function f(ua, Ma, e, l, m,
    ! k). This involves calculating the integral I = int(f(ua),
    ! ua->0,2pi), where ua is the true anomaly.  Along with ua, the
    ! corresponding mean anomaly Ma(ua) is passed into f
    !
    ! This implementation assumes that f(ua) is an even function,
    ! scaled to order-unity

    ! Set up initial values

    N = 1

    I = PI*(f(0._WP, 0._WP, e) + f(PI, PI, e))/2

    ! Refine the quadrature by repeatedly doubling N

    e_fac = sqrt((1._WP - e)*(1._WP + e))

    refine_loop : do

       ! Update the value of the quadrature

       S = 0._WP

       !$OMP PARALLEL DO REDUCTION(+:S) PRIVATE(ua,Ea,Ma) SCHEDULE(static, 4)
       update_loop : do j = 1, N

          ! Add contributions to the sum

          ua = PI*(2*j-1.5_WP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)  
          S = S + f(ua, Ma, e)
        
          ua = PI*(2*j-1.0_WP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)  
          S = S + f(ua, Ma, e)
        
          ua = PI*(2*j-0.5_WP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)  
          S = S + f(ua, Ma, e)
        
       end do update_loop

       dI = PI*S/(4*N) - 0.75_WP*I
       I = I + dI

       N = 4*N

       ! Check for convergence

       $ASSERT(N <= N_MAX,Too many iterations)

       dI_conv = MAX(tol, 2._WP*EPSILON(0._WP))

       if (abs(dI) < dI_conv) exit refine_loop

    end do refine_loop

    ! Double I (since we only integrate from 0 to pi)

    I = 2*I

    ! Finish

    return

  end function orbit_quad

end module gyre_tide_util
