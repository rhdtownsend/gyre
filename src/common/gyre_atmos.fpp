! Module   : gyre_atmos
! Purpose  : atmosphere utility routines
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

module gyre_atmos

  ! Uses

  use core_kinds

  use gyre_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  private

  public :: atmos_wavenumber
  public :: eval_atmos_cutoff_freqs
  public :: eval_atmos_coeffs_unno
  public :: eval_atmos_coeffs_jcd

  ! Procedures

contains

  function atmos_wavenumber (V_g, As, c_1, omega_c, l) result (lambda)

    real(WP)                :: V_g
    real(WP), intent(in)    :: As
    real(WP), intent(in)    :: c_1
    complex(WP), intent(in) :: omega_c
    integer, intent(in)     :: l
    complex(WP)             :: lambda

    real(WP)    :: omega_c_cutoff_lo
    real(WP)    :: omega_c_cutoff_hi
    complex(WP) :: gamma
    complex(WP) :: sgamma

    ! Calculate the radial wavenumber in the atmosphere

    if(AIMAG(omega_c) == 0._WP) then

       ! Calculate cutoff frequencies

       call eval_atmos_cutoff_freqs(V_g, As, c_1, l, omega_c_cutoff_lo, omega_c_cutoff_hi)

       ! Evaluate the wavenumber

       gamma = -4._WP*V_g*c_1*(omega_c**2 - omega_c_cutoff_lo**2)*(omega_c**2 - omega_c_cutoff_hi**2)/omega_c**2

       if(ABS(REAL(omega_c)) > omega_c_cutoff_hi) then

          ! Acoustic waves

          lambda = 0.5_WP*((V_g + As - 2._WP) - SQRT(gamma))

       elseif(ABS(REAL(omega_c)) < omega_c_cutoff_lo) then

          ! Gravity waves

          lambda = 0.5_WP*((V_g + As - 2._WP) + SQRT(gamma))

       else

          ! Evanescent

          lambda = 0.5_WP*((V_g + As - 2._WP) - SQRT(gamma))

       endif

    else

       ! Evaluate the wavenumber

       gamma = (As - V_g + 4._WP)**2 + 4*(l*(l+1)/(c_1*omega_c**2) - V_g)*(c_1*omega_c**2 - As)
       sgamma = SQRT(gamma)

       if(AIMAG(omega_c) > 0._WP) then

          ! Decaying oscillations; choose the wave with diverging
          ! energy density (see Townsend 2000b)

          if(REAL(sgamma) > 0._WP) then
             lambda = 0.5_WP*((V_g + As - 2._WP) + sgamma)
          else
             lambda = 0.5_WP*((V_g + As - 2._WP) - sgamma)
          endif

       else

          ! Growing oscillations; choose the wave with non-diverging
          ! energy density (see Townsend 2000b)

          if(REAL(sgamma) > 0._WP) then
             lambda = 0.5_WP*((V_g + As - 2._WP) - sgamma)
          else
             lambda = 0.5_WP*((V_g + As - 2._WP) + sgamma)
          endif

       endif

    end if

    ! Finish

    return

  end function atmos_wavenumber

!****

  subroutine eval_atmos_cutoff_freqs (V_g, As, c_1, l, omega_c_cutoff_lo, omega_c_cutoff_hi)

    real(WP), intent(in)  :: V_g
    real(WP), intent(in)  :: As
    real(WP), intent(in)  :: c_1
    integer, intent(in)   :: l
    real(WP), intent(out) :: omega_c_cutoff_lo
    real(WP), intent(out) :: omega_c_cutoff_hi

    real(WP) :: a
    real(WP) :: b
    real(WP) :: c

    ! Evaluate the atmospheric cutoff frequencies from the supplied coefficients

    a = -4._WP*V_g*c_1**2
    b = ((As - V_g + 4._WP)**2 + 4._WP*V_g*As + 4._WP*l*(l+1))*c_1
    c = -4._WP*l*(l+1)*As

    omega_c_cutoff_lo = SQRT((-b + SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    omega_c_cutoff_hi = SQRT((-b - SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    
    $ASSERT(omega_c_cutoff_hi >= omega_c_cutoff_lo,Incorrect cutoff frequency ordering)

    ! Finish

    return

  end subroutine eval_atmos_cutoff_freqs

!****
  
  subroutine eval_atmos_coeffs_unno (cf, x_o, V_g, As, c_1)

    class(model_t), intent(in) :: cf
    real(WP), intent(in)       :: x_o
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients (Unno et al. formulation)

    V_g = cf%V(x_o)/cf%Gamma_1(x_o)
    As = cf%As(x_o)
    c_1 = cf%c_1(x_o)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_unno

!****
  
  subroutine eval_atmos_coeffs_jcd (cf, x_o, V_g, As, c_1)

    class(model_t), intent(in) :: cf
    real(WP), intent(in)       :: x_o
    real(WP), intent(out)      :: V_g
    real(WP), intent(out)      :: As
    real(WP), intent(out)      :: c_1

    ! Evaluate atmosphere coefficients (JCD formulation)

    V_g = cf%V(x_o)/cf%Gamma_1(x_o)
    As = cf%V(x_o)*(1._WP-1._WP/cf%Gamma_1(x_o))
    c_1 = cf%c_1(x_o)

    ! Finish

    return

  end subroutine eval_atmos_coeffs_jcd

end module gyre_atmos
