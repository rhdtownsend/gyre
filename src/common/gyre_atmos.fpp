! Module   : gyre_atmos
! Purpose  : atmosphere utility routines
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

  ! Procedures

contains

  function atmos_beta_r_ (V_g, As, c_1, omega, lambda) result (beta)

    real(WP)             :: V_g
    real(WP), intent(in) :: As
    real(WP), intent(in) :: c_1
    real(WP), intent(in) :: omega
    real(WP), intent(in) :: lambda
    real(WP)             :: beta

    real(WP)      :: psi2
    logical, save :: warned = .FALSE.

    ! Calculate the atmospheric radial wavenumber, as defined by
    ! [Tow2000b] (real frequencies)

    psi2 = (As - V_g + 4._WP)**2 + 4._WP*(lambda/(c_1*omega**2) - V_g)*(c_1*omega**2 - As)

    if (psi2 < 0._WP) then

       if (.NOT. warned) then
          $WARN(WARNING: Discarding imaginary part of atmospheric radial wavenumber)
          warned = .TRUE.
       endif
       
       beta = 0.5_WP*(V_g + As - 2._WP)

    else

       beta = 0.5_WP*(V_g + As - 2._WP - SQRT(psi2))

    endif

    ! Finish

    return

  end function atmos_beta_r_

  !****

  function atmos_beta_c_ (V_g, As, c_1, omega, lambda) result (beta)

    real(WP)                :: V_g
    real(WP), intent(in)    :: As
    real(WP), intent(in)    :: c_1
    complex(WP), intent(in) :: omega
    complex(WP), intent(in) :: lambda
    complex(WP)             :: beta

    complex(WP) :: psi2
    complex(WP) :: psi

    ! Calculate the atmospheric radial wavenumber, as defined by
    ! [Tow2000b] (complex frequencies)

    psi2 = (As - V_g + 4._WP)**2 + 4._WP*(lambda/(c_1*omega**2) - V_g)*(c_1*omega**2 - As)
    psi = SQRT(psi2)

    ! Adjust the sign of psi to choose the correct solution

    if (AIMAG(omega) > 0._WP) then

       ! Decaying mode; choose the solution with diverging energy density

       if (REAL(psi) < 0._WP) psi = -psi

    elseif (AIMAG(omega) < 0._WP) then

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

          if (REAL(c_1*omega**2) < As) then
          
             ! Gravity wave; choose the solution with downward phase
             ! velocity

             if (AIMAG(psi) < 0._WP) psi = -psi

          else

             ! Acoustic wave; choose the solution with upward phase
             ! velocity

             if (AIMAG(psi) > 0._WP) psi = -psi

          end if

       end if

    end if

    beta = 0.5_WP*(V_g + As - 2._WP + psi)
             
    ! Finish

    return

  end function atmos_beta_c_

  !****

  subroutine eval_atmos_cutoff_freqs (V_g, As, c_1, lambda, omega_cutoff_lo, omega_cutoff_hi)

    real(WP), intent(in)  :: V_g
    real(WP), intent(in)  :: As
    real(WP), intent(in)  :: c_1
    real(WP), intent(in)  :: lambda
    real(WP), intent(out) :: omega_cutoff_lo
    real(WP), intent(out) :: omega_cutoff_hi

    real(WP) :: a
    real(WP) :: b
    real(WP) :: c

    ! Evaluate the atmospheric cutoff frequencies from the supplied coefficients

    a = -4._WP*V_g*c_1**2
    b = ((As - V_g + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda)*c_1
    c = -4._WP*lambda*As

    omega_cutoff_lo = SQRT((-b + SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    omega_cutoff_hi = SQRT((-b - SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    
    $ASSERT(omega_cutoff_hi >= omega_cutoff_lo,Incorrect cutoff frequency ordering)

    ! Finish

    return

  end subroutine eval_atmos_cutoff_freqs

  !****
  
  subroutine eval_atmos_coeffs_unno (ml, pt, V_g, As, c_1)

    class(model_t), intent(in) :: ml
    type(point_t), intent(in)  :: pt
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients ([Unn1989] formulation)

    V_g = ml%V_2(pt)*pt%x**2/ml%Gamma_1(pt)
    As = ml%As(pt)
    c_1 = ml%c_1(pt)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_unno

  !****
  
  subroutine eval_atmos_coeffs_jcd (ml, pt, V_g, As, c_1)

    class(model_t), intent(in) :: ml
    type(point_t), intent(in)  :: pt
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients ([Chr2008] formulation)

    V_g = ml%V_2(pt)*pt%x**2/ml%Gamma_1(pt)
    As = ml%V_2(pt)*pt%x**2*(1._WP-1._WP/ml%Gamma_1(pt))
    c_1 = ml%c_1(pt)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_jcd

end module gyre_atmos
