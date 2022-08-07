! Program  : gyre_tidal_coeff
! Purpose  : tide-related coefficient evaluation
!
! Copyright 2018-2022 Rich Townsend & The GYRE Team
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

module gyre_tidal_coeff

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_freq
  use gyre_func
  use gyre_grid
  use gyre_math
  use gyre_model
  use gyre_orbit_par
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  real(WP), parameter :: TOL = 1E-12_WP

  ! Access specifiers

  private

  public :: tidal_Phi_T
  public :: tidal_dPhi_T
  public :: tidal_Omega_orb
  public :: tidal_R_a
  public :: tidal_c
  public :: secular_G_1
  public :: secular_G_2
  public :: secular_G_3
  public :: secular_G_4
  public :: hansen_X
  public :: hansen_X_QP

  ! Procedures

contains

  function tidal_Phi_T (ml, or_p, x, l, m, k) result (Phi_T)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    real(WP), intent(in)          :: x
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: Phi_T

    real(WP) :: R_a
    real(WP) :: eps_T
    real(WP) :: c

    ! Evaluate the tidal forcing potential Phi_T at x, in units of G*M/R

    R_a = tidal_R_a(ml, or_p)
    eps_T = R_a**3*or_p%q

    c = tidal_c(ml, or_p, l, m, k)

    Phi_T = -eps_T/sqrt(4._WP*PI)*c*x**l

    ! Finish

    return

  end function tidal_Phi_T

  !****

  function tidal_dPhi_T (ml, or_p, x, l, m, k) result (dPhi_T)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    real(WP), intent(in)          :: x
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: dPhi_T

    real(WP) :: R_a
    real(WP) :: eps_T
    real(WP) :: c

    ! Evaluate the tidal forcing potential gradient dPhi_T/dx at x, in
    ! units of G*M/R

    R_a = tidal_R_a(ml, or_p)
    eps_T = R_a**3*or_p%q

    c = tidal_c(ml, or_p, l, m, k)

    dPhi_T = -eps_T/sqrt(4._WP*PI)*c*l*x**(l-1)

    ! Finish

    return

  end function tidal_dPhi_T

  !****

  function tidal_Omega_orb (ml, or_p) result (Omega_orb)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    real(WP)                      :: Omega_orb

    ! Evaluate the dimensionless orbital frequency

    Omega_orb = or_p%Omega_orb/freq_scale(or_p%Omega_orb_units, ml)

    ! Finish

    return

  end function tidal_Omega_orb

  !****

  function tidal_R_a (ml, or_p) result (R_a)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    real(WP)                      :: R_a

    real(WP) :: Omega_orb

    ! Evaluate the ratio of the stellar radius to the semi-major axis

    Omega_orb = tidal_Omega_orb(ml, or_p)

    R_a = (Omega_orb**2/(1._WP + or_p%q))**(1._WP/3._WP)

    ! Finish

    return

  end function tidal_R_a

  !****

  function tidal_c (ml, or_p, l, m, k) result (c)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: c

    real(WP) :: R_a
    real(WP) :: X

    ! Evaluate the tidal potential coefficient. This differs from the
    ! definition e.g. in eq. (5) of [Willems:2010aa], due to the
    ! different spherical harmonic definition

    R_a = tidal_R_a(ml, or_p)
    X = hansen_X(or_p, -(l+1), -m, -k)

    c = (4._WP*PI/(2*l+1))*R_a**(l-2)*REAL(CONJG(spherical_Y(l, m, HALFPI, 0._WP)))*X

    ! Finish

    return

  end function tidal_c

  !****

  function secular_G_1 (ml, or_p, l, m, k) result (G_1)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: G_1

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_1m1
    real(WP) :: X_1p1
    real(WP) :: X_2m1
    real(WP) :: X_2p1

    ! Evaluate the secular evolution coefficient

    c = tidal_c(ml, or_p, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_1m1 = hansen_X(or_p, -(l+1), -m-1, -k)
    X_1p1 = hansen_X(or_p, -(l+1), -m+1, -k)

    X_2m1 = hansen_X(or_p, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(or_p, -(l+2), -m+1, -k)

    associate (e => or_p%e)
      G_1 = c*Y* &
           (0.5_WP*(l+1)*(X_2m1 + X_2p1) + 0.5_WP*m*(X_2m1 - X_2p1) + &
           0.5_WP*m/(1._WP - e**2)*(X_1m1 - X_1p1))*sqrt(1._WP - e**2)/e
    end associate

    ! Finish

    return

  end function secular_G_1
  
  !****

  function secular_G_2 (ml, or_p, l, m, k) result (G_2)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: G_2

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_2m1
    real(WP) :: X_2p1
    real(WP) :: X_3

    ! Evaluate the secular evolution coefficient

    c = tidal_c(ml, or_p, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_2m1 = hansen_X(or_p, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(or_p, -(l+2), -m+1, -k)

    X_3 = hansen_X(or_p, -(l+3), -m, -k)

    associate (e => or_p%e)
      G_2 = -2._WP*c*Y* &
           (0.5_WP*(l+1)*e*(X_2m1 - X_2p1) + m*(1._WP - e**2)*X_3)/sqrt(1._WP - e**2)
    end associate

    ! Finish

    return

  end function secular_G_2

  !****

  function secular_G_3 (ml, or_p, l, m, k) result (G_3)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: G_3

    real(WP) :: c
    real(WP) :: Y
    real(WP) :: X_1
    real(WP) :: X_2m1
    real(WP) :: X_2p1
    real(WP) :: X_3

    ! Evaluate the secular evolution coefficient

    c = tidal_c(ml, or_p, l, m, k)
    Y = REAL(spherical_Y(l, m, HALFPI, 0._WP))

    X_1 = hansen_X(or_p, -(l+1), -m, -k)

    X_2m1 = hansen_X(or_p, -(l+2), -m-1, -k)
    X_2p1 = hansen_X(or_p, -(l+2), -m+1, -k)

    X_3 = hansen_X(or_p, -(l+3), -m, -k)

    associate (e => or_p%e)
      G_3 = -c*Y* &
           (0.5_WP*(l+1)*e*(X_2m1 - X_2p1) + m*(1._WP - e**2)*X_3 - m*X_1)*sqrt(1._WP - e**2)/e
    end associate

    ! Finish

    return

  end function secular_G_3

  !****

  function secular_G_4 (ml, or_p, l, m, k) result (G_4)

    class(model_t), intent(in)    :: ml
    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: l
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: G_4

    real(WP) :: R_a
    real(WP) :: c

    ! Evaluate the secular evolution coefficient

    R_a = tidal_R_a(ml, or_p)
    c = tidal_c(ml, or_p, l, m, k)

    G_4 = m*(2*l+1)/(4*PI)*(R_a)**(-l+2)*c**2

    ! Finish

    return

  end function secular_G_4

  !****

  function hansen_X (or_p, n, m, k) result (X)

    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: n
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: X

    real(WP) :: A
    real(WP) :: S

    ! Evaluate the Hansen coefficient X_nmk

    associate (e => or_p%e)

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

    end associate

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

  function hansen_X_QP (or_p, n, m, k) result (X)

    type(orbit_par_t), intent(in) :: or_p
    integer, intent(in)           :: n
    integer, intent(in)           :: m
    integer, intent(in)           :: k
    real(WP)                      :: X

    real(WP) :: A
    real(WP) :: S

    ! Evaluate the Hansen coefficient X_nmk

    associate (e => or_p%e)

      if (e**2 > EPSILON(0._WP)) then

         A = (1._WP - e**2)**(n+1.5_WP)/TWOPI
         S = MAX((1._WP - e)**(-n-2), (1._WP + e)**(-n-2))

         X = REAL(A*S*orbit_quad_QP(hansen_X_f_, REAL(e, QP), REAL(TOL, QP)), WP)

      else

         if (m == k) then
            X = 1._WP
         else
            X = 0._WP
         endif

      endif

    end associate

    ! Finish

    return

  contains

    function hansen_X_f_ (ua, Ma, e) result (f)

      real(QP), intent(in) :: ua
      real(QP), intent(in) :: Ma
      real(QP), intent(in) :: e
      real(QP)             :: f

      ! Set up the integrand

      f = cos(m*ua - k*Ma)/(S*(1._QP + e*cos(ua))**(n+2))

      ! Finish

      return

    end function hansen_X_f_

  end function hansen_X_QP

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

       ! Parallelization disabled for now to give reproducible results
       ! !$OMP PARALLEL DO REDUCTION(+:S) PRIVATE(ua,Ea,Ma) SCHEDULE(static, 4)
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

  !****

  function orbit_quad_QP (f, e, tol) result (I)

    interface
       function f (ua, Ma, e)
         use core_kinds
         real(QP), intent(in) :: ua
         real(QP), intent(in) :: Ma
         real(QP), intent(in) :: e
         real(QP)             :: f
       end function f
    end interface
    real(QP), intent(in) :: e
    real(QP), intent(in) :: tol
    real(QP)             :: I

    integer, parameter :: N_MAX = 2**24
    real(QP), parameter :: PI = ACOS(-1._QP)

    integer  :: N
    real(QP) :: e_fac
    real(QP) :: S
    integer  :: j
    real(QP) :: ua
    real(QP) :: Ea
    real(QP) :: Ma
    real(QP) :: dI
    real(QP) :: dI_conv

    ! Perform "orbital quadrature" on the function f(ua, Ma, e, l, m,
    ! k). This involves calculating the integral I = int(f(ua),
    ! ua->0,2pi), where ua is the true anomaly.  Along with ua, the
    ! corresponding mean anomaly Ma(ua) is passed into f
    !
    ! This implementation assumes that f(ua) is an even function,
    ! scaled to order-unity

    ! Set up initial values

    N = 1

    I = PI*(f(0._QP, 0._QP, e) + f(PI, PI, e))/2

    ! Refine the quadrature by repeatedly doubling N

    e_fac = sqrt((1._QP - e)*(1._QP + e))

    refine_loop : do

       ! Update the value of the quadrature

       S = 0._QP

       !$OMP PARALLEL DO REDUCTION(+:S) PRIVATE(ua,Ea,Ma) SCHEDULE(static, 4)
       update_loop : do j = 1, N

          ! Add contributions to the sum

          ua = PI*(2*j-1.5_QP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)
          S = S + f(ua, Ma, e)

          ua = PI*(2*j-1.0_QP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)
          S = S + f(ua, Ma, e)

          ua = PI*(2*j-0.5_QP)/(2*N)
          Ea = atan2(e_fac*sin(ua), e + cos(ua))
          Ma = Ea - e*sin(Ea)
          S = S + f(ua, Ma, e)

       end do update_loop

       dI = PI*S/(4*N) - 0.75_QP*I
       I = I + dI

       N = 4*N

       ! Check for convergence

       $ASSERT(N <= N_MAX,Too many iterations)

       dI_conv = MAX(tol, 2._QP*EPSILON(0._QP))

       if (abs(dI) < dI_conv) exit refine_loop

    end do refine_loop

    ! Double I (since we only integrate from 0 to pi)

    I = 2*I

    ! Finish

    return

  end function orbit_quad_QP

end module gyre_tidal_coeff
