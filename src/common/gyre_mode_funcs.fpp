! Module   : gyre_mode_funcs
! Purpose  : support functions for gyre_mode
!
! Copyright 2013 Rich Townsend
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
  use gyre_modepar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: xi_r
  public :: xi_h
  public :: eul_phi
  public :: deul_phi_dx
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
  public :: F_j_rey
  public :: Yt_1
  public :: Yt_2
  public :: I_0
  public :: I_1

  ! Procedures

contains

  function xi_r (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: xi_r

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the radial displacement perturbation at x, in units of
    ! R_star

    associate (l => mp%l)

      if (l /= 1) then

         if (x /= 0._WP) then
            xi_r = y(1)*x**(l-1)
         else
            xi_r = 0._WP
         endif

      else

         xi_r = y(1)

      endif

    end associate

    ! Finish

    return

  end function xi_r

!****

  function xi_h (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: xi_h

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the horizontal displacement perturbation at x, in
    ! units of R_star

    associate (c_1 => ml%c_1(x), l => mp%l, omega_c => ml%omega_c(x, mp%m, omega))

      if (l /= 0) then

         if(l /= 1) then

            if (x /= 0._WP) then
               xi_h = y(2)*x**(l-1)/(c_1*omega_c**2)
            else
               xi_h = 0._WP
            end if

         else
            
            xi_h = y(2)/(c_1*omega_c**2)

         endif

      else

         xi_h = 0._WP

      end if

    end associate

    ! Finish

    return

  end function xi_h

!****

  function eul_phi (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: eul_phi

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the Eulerian gravitational potential perturbation at
    ! x, in units of G M_star / R_star

    associate (c_1 => ml%c_1(x), l => mp%l)

      eul_phi = y(3)*x**l/c_1

    end associate

    ! Finish

    return

  end function eul_phi

!****

  function deul_phi (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: deul_phi_dx

    $CHECK_BOUNDS(SIZE(y),6)
    
    ! Calculate the Eulerian potential gradient (gravity) perturbation
    ! at x, in units of G M_star / R_star**2

    associate (c_1 => ml%c_1(x), l => mp%l)

      if (l /= 1) then

         if (x /= 0._WP) then
            deul_phi_dx = y(4)*x**(l-1)/c_1
         else
            deul_phi_dx = 0._WP
         end if

      else

         deul_phi_dx = y(4)/c_1

      end if

    end associate

    ! Finish

    return

  end function deul_phi

!****

  function lag_S (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_S

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian specific entropy perturbation at x, in units
    ! of c_p

    associate (l => mp%l)

      if (x /= 0._WP) then
         lag_S = y(5)*x**(l-2)
      else
         lag_S = 0._WP
      endif

    end associate
         
    ! Finish

    return

  end function lag_S

!****

  function lag_L (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_L

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian luminosity perturbation at x, in units
    ! of L_star

    associate (l => mp%l)

      lag_L = y(6)*x**(l+1)

    end associate

    ! Finish

    return

  end function lag_L

!****

  function eul_P (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: eul_P

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Eulerian pressure perturbation at x, in units of
    ! P

    associate (V => ml%V(x))

      if (x /= 0._WP) then
         eul_P = lag_P(ml, mp, omega, x, y) + V*xi_r(ml, mp, omega, x, y)/x
      else
         eul_P = delp(ml, mp, omega, x, y)
      endif

    end associate

    ! Finish

    return

  end function eul_P

!****

  function lag_P (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_P

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian pressure perturbation at x, in units of
    ! P

    associate (V => ml%V(x), pi_c => ml%pi_c(), l => mp%l)

      if(l > 0) then

         if (x /= 0._WP) then
            lag_P = V*(y(2) - y(1) - y(3))*x**(l-2)
         else
            lag_P = 0._WP
         end if

      else

         if (x /= 0._WP) then
            lag_P = V*(y(2) - y(1) - y(3))*x**(l-2)
         else
            lag_P = pi_c*(y(2) - y(1) - y(3))
         endif

      endif

    end associate

    ! Finish

    return

  end function lag_P

!****

  function eul_rho (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: eul_rho

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Eulerian density perturbation at x, in units of
    ! rho

    associate (V_g => ml%V(x)/ml%Gamma_1(x), As => ml%As(x))

      if (x /= 0._WP) then
         eul_rho = lag_rho(ml, mp, omega, x, y) + (V_g + As)*xi_r(ml, mp, omega, x, y)/x
      else
         eul_rho = lag_rho(ml, mp, omega, x, y)
      endif

    end associate

    ! Finish

    return

  end function eul_rho

!****

  function lag_rho (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_rho

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian density perturbation at x, in units of
    ! rho. This expression implements eqn. 13.83 of [Unn1989]

    associate (Gamma_1 => ml%Gamma_1(x), delta => ml%delta(x))

      lag_rho = lag_P(ml, mp, omega, x, y)/Gamma_1 - delta*lag_S(ml, mp, omega, x, y)

    end associate

    ! Finish

    return

  end function lag_rho

!****

  function eul_T (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: eul_T

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian temperature perturbation at x, in units
    ! of T

    associate (V => ml%V(x), nabla => ml%nabla(x), nabla_ad => ml%nabla_ad(x))
      
      if (x /= 0._WP) then
         eul_T = lag_T(ml, mp, omega, x, y) + nabla*V*xi_r(ml, mp, omega, x, y)/x
      else
         eul_T = lag_T(ml, mp, omega, x, y)
      endif

    end associate

    ! Finish

    return

  end function eul_T

!****

  function lag_T (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: lag_T

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Lagrangian temperature perturbation at x, in units
    ! of T. This expression implements eqn. 13.84 of [Unn1989]

    associate (nabla_ad => ml%nabla_ad(x))
      
      lag_T = nabla_ad*lag_P(ml, mp, omega, x, y) + lag_S(ml, mp, omega, x, y)

    end associate

    ! Finish

    return

  end function lag_T

!****

  function dE_dx (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    real(WP)                    :: dE_dx

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the differential mode inertia at x, in units of M_star
    ! R_star**2. This expression is based on eqn. 3.139 of [Aer2010]

    associate(xi_r => xi_r(ml, mp, omega, x, y), xi_h => xi_h(ml, mp, omega, x, y), &
              U => ml%U(x), c_1 => ml%c_1(x), l => mp%l)
      dE_dx = (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*x**2/c_1
    end associate

    ! Finish

    return

  end function dE_dx

!****

  function dW_dx (ml, mp, omega, x, y)

    use gyre_evol_model

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    real(WP)                    :: dW_dx

    real(WP) :: t_dyn
    real(WP) :: t_kh

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

    associate(c_thm => ml%c_thm(x))

      dW_dx = -PI*AIMAG(CONJG(lag_T(ml, mp, omega, x, y))*lag_S(ml, mp, omega, x, y))*c_thm*x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_dx

!****

  function F_j_rey (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    real(WP)                    :: F_j_rey
    
    ! Calculate the angle-averaged angular momentum flux due to Reynolds stress, in units of G
    ! M_star**2/R_star**3.  This expression is based on eqn. 21 of [LeeSai1993]

    associate (c_1 => ml%c_1(x), U => ml%U(x), m => mp%m, omega_c => ml%omega_c(x, mp%m, omega))

      F_j_rey = -ABS(omega_c**2)*x*U*AIMAG(CONJG(xi_r(ml, mp, omega, x, y))*m*xi_h(ml, mp, omega, x, y))/(32._WP*PI**2*c_1)

    end associate

  end function F_j_rey

!****

  function Yt_1 (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: Yt_1

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Takata Y_1 function at x. This expression is based
    ! on eqn. 69 of [Tak2006b]

    associate (J => 1._WP-ml%U(x)/3._WP)

      Yt_1 = J*y(1) + (y(3) - y(4))/3._WP
      
    end associate

    ! Finish

    return

  end function Yt_1

!****

  function Yt_2 (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: Yt_2

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the Takata Y_2 function at x. This expression is based
    ! on eqn. 70 of [Tak2006b], divided through by V

    Yt_2 = y(2) - y(1) - y(3)

    ! Finish

    return

  end function Yt_2

!****

  function I_0 (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: I_0

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the I_0 integral at x, which should be zero for radial
    ! modes. This expression is based on eqn. 42 of [Tak2006a]

    associate(U => ml%U(x), c_1 => ml%c_1(x), l => mp%l)

      I_0 = x**(l+1)*(U*y(1) + y(4))/c_1

    end associate

    ! Finish

    return

  end function I_0

!****

  function I_1 (ml, mp, omega, x, y)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x
    complex(WP), intent(in)     :: y(:)
    complex(WP)                 :: I_1

    $CHECK_BOUNDS(SIZE(y),6)

    ! Calculate the I_0 integral at x, which should be zero for dipole
    ! modes. This expression is based on eqn. 43 of [Tak2006a]

    associate(U => ml%U(x), c_1 => ml%c_1(x), l => mp%l, &
              omega_c => ml%omega_c(x, mp%m, omega))

      I_1 = x**(l+2)*(c_1*omega_c**2*U*y(1) - U*y(2) + &
                  (U - c_1*omega_c**2 - 2._WP)*y(3) + (c_1*omega_c**2 - 1._WP)*y(4))/c_1**2

    end associate

    ! Finish

    return

  end function I_1

end module gyre_mode_funcs
