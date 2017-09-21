! Module   : gyre_atmos
! Purpose  : atmosphere utility routines
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_atmos

  ! Uses

  use core_kinds

  use gyre_point
  use gyre_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface atmos_beta
     module procedure atmos_beta_r_
     module procedure atmos_beta_c_
  end interface atmos_beta

  ! Access specifiers

  private

  public :: atmos_beta
  public :: eval_atmos_cutoff_freqs
  public :: eval_atmos_coeffs_unno
  public :: eval_atmos_coeffs_jcd
  public :: eval_atmos_coeffs_luan

  ! Procedures

contains

  function atmos_beta_r_ (V_g, As, U, c_1, omega, lambda) result (beta)

    real(WP)             :: V_g
    real(WP), intent(in) :: As
    real(WP), intent(in) :: U
    real(WP), intent(in) :: c_1
    real(WP), intent(in) :: omega
    real(WP), intent(in) :: lambda
    real(WP)             :: beta

    real(WP)      :: b_11
    real(WP)      :: b_12
    real(WP)      :: b_21
    real(WP)      :: b_22
    real(WP)      :: psi2
    logical, save :: warned = .FALSE.

    ! Calculate the atmospheric radial wavenumber, as defined by
    ! [Tow2000b] (real frequencies)

    b_11 = V_g - 3._WP
    b_12 = lambda/(c_1*omega**2) - V_g
    b_21 = c_1*omega**2 - As
    b_22 = As - U + 1._WP

    psi2 = (b_22 - b_11)**2 + 4._WP*b_12*b_21

    if (psi2 < 0._WP) then

       if (.NOT. warned) then
          $WARN(WARNING: Discarding imaginary part of atmospheric radial wavenumber)
          warned = .TRUE.
       endif
       
       beta = 0.5_WP*(b_11 + b_22)

    else

       beta = 0.5_WP*(b_11 + b_22 - SQRT(psi2))

    endif

    ! Finish

    return

  end function atmos_beta_r_

  !****

  function atmos_beta_c_ (V_g, As, U, c_1, omega, lambda) result (beta)

    real(WP)                :: V_g
    real(WP), intent(in)    :: As
    real(WP), intent(in)    :: U
    real(WP), intent(in)    :: c_1
    complex(WP), intent(in) :: omega
    complex(WP), intent(in) :: lambda
    complex(WP)             :: beta

    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: psi2
    complex(WP) :: psi

    ! Calculate the atmospheric radial wavenumber, as defined by
    ! [Tow2000b] (complex frequencies)

    b_11 = V_g - 3._WP
    b_12 = lambda/(c_1*omega**2) - V_g
    b_21 = c_1*omega**2 - As
    b_22 = As - U + 1

    psi2 = (b_22 - b_11)**2 + 4._WP*b_12*b_21
    psi = SQRT(psi2)

    ! Adjust the sign of psi to choose the correct solution

    if (AIMAG(omega) < 0._WP) then

       ! Decaying mode; choose the solution with diverging energy density

       if (REAL(psi) < 0._WP) psi = -psi

    elseif (AIMAG(omega) > 0._WP) then

       ! Growing mode; choose the solution with non-diverging energy
       ! density

       if (REAL(psi) > 0._WP) psi = -psi

    else

       ! Steady-state mode;

       $ASSERT(AIMAG(psi2) == 0._WP,Steady-state mode with complex propagation discriminant)

       if (REAL(psi2) > 0._WP) then

          ! Evanescent wave; choose the solution with non-diverging
          ! energy density

          psi = -psi

       else

          if (REAL(b_21) < 0._WP) then
          
             ! Gravity wave (c_1*omega**2 < As); choose the solution
             ! with downward phase velocity

             if (AIMAG(psi) < 0._WP) psi = -psi

          else

             ! Acoustic wave; choose the solution with upward phase
             ! velocity

             if (AIMAG(psi) > 0._WP) psi = -psi

          end if

       end if

    end if

    beta = 0.5_WP*(b_11 + b_22 + psi)
             
    ! Finish

    return

  end function atmos_beta_c_

  !****

  subroutine eval_atmos_cutoff_freqs (V_g, As, U, c_1, lambda, omega_cutoff_lo, omega_cutoff_hi)

    real(WP), intent(in)  :: V_g
    real(WP), intent(in)  :: As
    real(WP), intent(in)  :: U
    real(WP), intent(in)  :: c_1
    real(WP), intent(in)  :: lambda
    real(WP), intent(out) :: omega_cutoff_lo
    real(WP), intent(out) :: omega_cutoff_hi

    real(WP) :: a
    real(WP) :: b
    real(WP) :: c

    ! Evaluate the atmospheric cutoff frequencies from the supplied coefficients

    a = -4._WP*V_g*c_1**2
    b = ((As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda)*c_1
    c = -4._WP*lambda*As

    omega_cutoff_lo = SQRT((-b + SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    omega_cutoff_hi = SQRT((-b - SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    
    $ASSERT(omega_cutoff_hi >= omega_cutoff_lo,Incorrect cutoff frequency ordering)

    ! Finish

    return

  end subroutine eval_atmos_cutoff_freqs

  !****
  
  subroutine eval_atmos_coeffs_unno (ml, pt, V_g, As, U, c_1)

    class(model_t), intent(in) :: ml
    type(point_t), intent(in)  :: pt
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: U
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients ([Unn1989] formulation)

    V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
    As = ml%coeff(I_AS, pt)
    U = 0._WP
    c_1 = ml%coeff(I_C_1, pt)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_unno

  !****
  
  subroutine eval_atmos_coeffs_jcd (ml, pt, V_g, As, U, c_1)

    class(model_t), intent(in) :: ml
    type(point_t), intent(in)  :: pt
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: U
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients ([Chr2008] formulation)

    V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
    As = ml%coeff(I_V_2, pt)*pt%x**2*(1._WP-1._WP/ml%coeff(I_GAMMA_1, pt))
    U = 0._WP
    c_1 = ml%coeff(I_C_1, pt)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_jcd

  !****
  
  subroutine eval_atmos_coeffs_luan (ml, pt, V_g, As, U, c_1)

    class(model_t), intent(in) :: ml
    type(point_t), intent(in)  :: pt
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: U
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients (Luan formulation)

    V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
    As = 0._WP
    U = ml%coeff(I_U, pt)
    c_1 = ml%coeff(I_C_1, pt)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_luan

end module gyre_atmos
