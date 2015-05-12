! Module   : gyre_mode_funcs
! Purpose  : support functions for gyre_mode
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_mode_funcs

  ! Uses

  use core_kinds
  use gyre_constants

  use gyre_model
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: xi_r
  public :: xi_h
  public :: eul_phi
  public :: deul_phi
  public :: lag_S
  public :: lag_L
  public :: eul_P
  public :: lag_P
  public :: eul_rho
  public :: lag_rho
  public :: eul_T
  public :: lag_T
  public :: dE_dx
  public :: dW_dx
  public :: F_j
  public :: div_F_j
  public :: Yt_1
  public :: Yt_2
  public :: I_0
  public :: I_1

  ! Procedures

contains

  function xi_r (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: xi_r

    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the radial displacement perturbation at x, in units of
    ! R_star

    l_0 = rt%l_0(omega)

    if (l_0 /= 1._WP) then

       if (x /= 0._WP) then
          xi_r = y(1)*x**(l_0-1._WP)
       else
          xi_r = 0._WP
       endif

    else

       xi_r = y(1)

    endif

    ! Finish

    return

  end function xi_r

!****

  function xi_h (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: xi_h

    real(WP)    :: c_1
    complex(WP) :: l_0
    complex(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the horizontal displacement perturbation at x, in
    ! units of R_star

    c_1 = ml%c_1(x)

    l_0 = rt%l_0(omega)

    omega_c = rt%omega_c(x, omega)

    if (l_0 /= 0._WP) then

       if(l_0 /= 1._WP) then

          if (x /= 0._WP) then
             xi_h = y(2)*x**(l_0-1._WP)/(c_1*omega_c**2)
          else
             xi_h = 0._WP
          end if

       else
            
          xi_h = y(2)/(c_1*omega_c**2)

       endif

    else

       xi_h = 0._WP

    end if

    ! Finish

    return

  end function xi_h

!****

  function eul_phi (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: eul_phi

    real(WP)    :: c_1
    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the Eulerian gravitational potential perturbation at
    ! x, in units of G M_star / R_star

    c_1 = ml%c_1(x)

    l_0 = rt%l_0(omega)

    if (l_0 /= 0._WP) then

       if (x /= 0._WP) then
          eul_phi = y(3)*x**l_0/c_1
       else
          eul_phi = 0._WP
       endif

    else

       eul_phi = y(3)/c_1

    endif

    ! Finish

    return

  end function eul_phi

!****

  function deul_phi (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: deul_phi

    real(WP)    :: c_1
    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the Eulerian potential gradient (gravity) perturbation
    ! at x, in units of G M_star / R_star**2

    c_1 = ml%c_1(x)

    l_0 = rt%l_0(omega)

    if (l_0 /= 1._WP) then

       if (x /= 0._WP) then
          deul_phi = y(4)*x**(l_0-1._WP)/c_1
       else
          deul_phi = 0._WP
       end if
       
    else

       deul_phi = y(4)/c_1

    end if

    ! Finish

    return

  end function deul_phi

!****

  function lag_S (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: lag_S

    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian specific entropy perturbation at x, in
    ! units of c_p

    l_0 = rt%l_0(omega)

    if (x /= 0._WP) then
       lag_S = y(5)*x**(l_0-2._WP)
    else
       lag_S = 0._WP
    endif

    ! Finish

    return

  end function lag_S

!****

  function lag_L (ml, rt, omega, x, y)

    class(model_t), intent(in)  :: ml
    class(c_rot_t), intent(in)  :: rt
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_L

    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian luminosity perturbation at x, in units
    ! of L_star

    l_0 = rt%l_0(omega)

    if (x /= 0._WP) then
       lag_L = y(6)*x**(l_0+1._WP)
    else
       lag_L = 0._WP
    endif

    ! Finish

    return

  end function lag_L

!****

  function eul_P (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: eul_P

    real(WP) :: V_2

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Eulerian pressure perturbation at x, in units of
    ! P

    V_2 = ml%V_2(x)

    eul_P = lag_P(ml, rt, omega, x, y) + V_2*x*xi_r(ml, rt, omega, x, y)

    ! Finish

    return

  end function eul_P

!****

  function lag_P (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: lag_P

    real(WP)    :: V_2
    complex(WP) :: l_0

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian pressure perturbation at x, in units of
    ! P

    V_2 = ml%V_2(x)

    l_0 = rt%l_0(omega)

    if (l_0 /= 0._WP) then

       if (x /= 0._WP) then
          lag_P = V_2*(y(2) - y(1) - y(3))*x**l_0
       else
          lag_P = 0._WP
       end if

    else

       lag_P = V_2*(y(2) - y(1) - y(3))

    endif

    ! Finish

    return

  end function lag_P

!****

  function eul_rho (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: eul_rho

    real(WP) :: D

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Eulerian density perturbation at x, in units of
    ! rho

    D = ml%D(x)

    if (x /= 0._WP) then
       eul_rho = lag_rho(ml, rt, omega, x, y) - D*xi_r(ml, rt, omega, x, y)/x
    else
       eul_rho = lag_rho(ml, rt, omega, x, y)
    endif

    ! Finish

    return

  end function eul_rho

!****

  function lag_rho (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: lag_rho

    real(WP) :: Gamma_1
    real(WP) :: delta

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian density perturbation at x, in units of
    ! rho. This expression implements eqn. 13.83 of [Unn1989]

    Gamma_1 = ml%Gamma_1(x)
    delta = ml%delta(x)

    lag_rho = lag_P(ml, rt, omega, x, y)/Gamma_1 - delta*lag_S(ml, rt, omega, x, y)

    ! Finish

    return

  end function lag_rho

!****

  function eul_T (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: eul_T

    real(WP) :: V_2
    real(WP) :: nabla

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian temperature perturbation at x, in units
    ! of T

    V_2 = ml%V_2(x)
    nabla = ml%nabla(x)
      
    eul_T = lag_T(ml, rt, omega, x, y) + nabla*V_2*x*xi_r(ml, rt, omega, x, y)

    ! Finish

    return

  end function eul_T

!****

  function lag_T (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: lag_T

    real(WP) :: nabla_ad

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian temperature perturbation at x, in units
    ! of T. This expression implements eqn. 13.84 of [Unn1989]

    nabla_ad = ml%nabla_ad(x)
      
    lag_T = nabla_ad*lag_P(ml, rt, omega, x, y) + lag_S(ml, rt, omega, x, y)

    ! Finish

    return

  end function lag_T

!****

  function dE_dx (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    real(WP)                   :: dE_dx

    complex(WP) :: xi_r_
    complex(WP) :: xi_h_
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: lambda

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the differential mode inertia at x, in units of M_star
    ! R_star**2. This expression is based on eqn. 3.139 of [Aer2010]

    xi_r_ = xi_r(ml, rt, omega, x, y)
    xi_h_ = xi_h(ml, rt, omega, x, y)

    U = ml%U(x)
    c_1 = ml%c_1(x)

    lambda = rt%lambda(x, omega)

    dE_dx = (ABS(xi_r_)**2 + ABS(lambda)*ABS(xi_h_)**2)*U*x**2/c_1

    ! Finish

    return

  end function dE_dx

!****

  function dW_dx (ml, rt, omega, x, y)

    use gyre_evol_model

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    real(WP)                   :: dW_dx

    real(WP) :: t_dyn
    real(WP) :: t_kh
    real(WP) :: c_thm

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the differential work at x, in units of G
    ! M_star**2/R_star.  This expression is based on eqn. 25.9 of
    ! [Unn1989]

    select type (ml)
    class is (evol_model_t)
       t_dyn = SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    c_thm = ml%c_thm(x)

    dW_dx = -PI*AIMAG(CONJG(lag_T(ml, rt, omega, x, y))*lag_S(ml, rt, omega, x, y))*c_thm*x**2*t_dyn/t_kh

    ! Finish

    return

  end function dW_dx

!****

  function F_j (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    real(WP)                   :: F_j

    real(WP)    :: c_1
    real(WP)    :: U
    integer     :: m
    complex(WP) :: omega_c
    
    ! Calculate the angle-averaged angular momentum flux due to
    ! Reynolds stress, in units of G M_star**2/R_star**3.  This
    ! expression is based on eqn. 21 of [LeeSai1993]

    c_1 = ml%c_1(x)
    U = ml%U(x)

    m = rt%mp%m

    omega_c = rt%omega_c(x, omega)

    F_j = ABS(omega_c**2)*x*U*AIMAG(CONJG(xi_r(ml, rt, omega, x, y))*m*xi_h(ml, rt, omega, x, y))/(32._WP*PI**2*c_1)

    ! Finish

    return

  end function F_j

!****

  function div_F_j (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    real(WP)                   :: div_F_j

    real(WP)    :: V_2
    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: U
    real(WP)    :: c_1
    integer     :: m
    complex(WP) :: omega_c
    real(WP)    :: div_F_j_wave
    real(WP)    :: div_F_j_NA
    
    ! Calculate the divergence of the angle-averaged angular momentum
    ! flux due to Reynolds stress, in units of G M_star**2/R_star**4.
    ! This expression is based on eqns. 8-10 of [And1983]

    V_2 = ml%V_2(x)
    V_g = V_2*x**2/ml%Gamma_1(x)
    As = ml%As(x)
    U = ml%U(x)
    c_1 = ml%c_1(x)
    
    m = rt%mp%m

    omega_c = rt%omega_c(x, omega)

    div_F_j_wave = REAL(omega_c)*AIMAG(omega_c)*m*(U/c_1)*(V_g*c_1*REAL(omega_c)**2*ABS(xi_h(ml, rt, omega, x, y))**2 + &
                                                           As/(c_1*REAL(omega_c)**2)*ABS(xi_r(ml, rt, omega, x, y))**2)/(16._WP*PI**2)
    div_F_j_NA = m*AIMAG(CONJG(lag_P(ml, rt, omega, x, y))*lag_rho(ml, rt, omega, x, y))*U/(32._WP*PI**2*c_1**2*V_2)

    div_F_j = div_F_j_wave + div_F_j_NA

    ! Finish

    return

  end function div_F_j

!****

  function Yt_1 (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: Yt_1

    real(WP) :: J

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Takata Y_1 function at x. This expression is based
    ! on eqn. 69 of [Tak2006b]

    J = 1._WP - ml%U(x)/3._WP

    Yt_1 = J*y(1) + (y(3) - y(4))/3._WP
      
    ! Finish

    return

  end function Yt_1

!****

  function Yt_2 (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: Yt_2

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Takata Y_2 function at x. This expression is based
    ! on eqn. 70 of [Tak2006b], divided through by V

    Yt_2 = y(2) - y(1) - y(3)

    ! Finish

    return

  end function Yt_2

!****

  function I_0 (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: I_0

    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: l_0
    
    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the I_0 integral at x, which should be zero for radial
    ! modes. This expression is based on eqn. 42 of [Tak2006a]

    U = ml%U(x)
    c_1 = ml%c_1(x)

    l_0 = rt%l_0(omega)

    if (x /= 0._WP) then
       I_0 = x**(l_0+1._WP)*(U*y(1) + y(4))/c_1
    else
       I_0 = 0._WP
    endif

    ! Finish

    return

  end function I_0

!****

  function I_1 (ml, rt, omega, x, y)

    class(model_t), intent(in) :: ml
    class(c_rot_t), intent(in) :: rt
    complex(WP), intent(in)    :: omega
    real(WP), intent(in)       :: x
    complex(WP), intent(in)    :: y(:)
    complex(WP)                :: I_1

    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: l_0
    complex(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the I_0 integral at x, which should be zero for dipole
    ! modes. This expression is based on eqn. 43 of [Tak2006a]

    U = ml%U(x)
    c_1 = ml%c_1(x)

    l_0 = rt%l_0(omega)

    omega_c = rt%omega_c(x, omega)

    if (x /= 0._WP) then
       I_1 = x**(l_0+2._WP)*(c_1*omega_c**2*U*y(1) - U*y(2) + &
            (U - c_1*omega_c**2 - 2._WP)*y(3) + (c_1*omega_c**2 - 1._WP)*y(4))/c_1**2
    else
       I_1 = 0._WP
    endif
      
    ! Finish

    return

  end function I_1

end module gyre_mode_funcs
