! Program  : gyre_units
! Purpose  : unit conversion routines
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

module gyre_units

  ! Uses

  use core_kinds
  use core_constants

  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_poly_base_coeffs
  use gyre_hom_base_coeffs
  use gyre_oscpar
  use gyre_ad_bound

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: freq_scale

contains

  function freq_scale (bc, op, x_o, freq_units)

    class(base_coeffs_t), intent(in) :: bc
    type(oscpar_t), intent(in)       :: op
    real(WP), intent(in)             :: x_o
    character(LEN=*), intent(in)     :: freq_units
    real(WP)                         :: freq_scale

    ! Calculate the scale factor to convert a dimensionless angular
    ! frequency to a dimensioned frequency

    select type (bc)
    class is (evol_base_coeffs_t)
       freq_scale = evol_freq_scale(bc, op, x_o, freq_units)
    class is (poly_base_coeffs_t)
       freq_scale = poly_freq_scale(freq_units)
    class is (hom_base_coeffs_t)
       freq_scale = hom_freq_scale(freq_units)
    class default
       $ABORT(Invalid bc type)
    end select

    ! Finish

    return

  contains

    function evol_freq_scale (bc, op, x_o, freq_units) result (freq_scale)

      class(evol_base_coeffs_t), intent(in) :: bc
      type(oscpar_t), intent(in)            :: op
      real(WP), intent(in)                  :: x_o
      character(LEN=*), intent(in)          :: freq_units
      real(WP)                              :: freq_scale

      real(WP) :: omega_cutoff_lo
      real(WP) :: omega_cutoff_hi

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case(freq_units)
      case('NONE')
         freq_scale = 1._WP
      case('HZ')
         freq_scale = 1._WP/(TWOPI*SQRT(bc%R_star**3/(bc%G*bc%M_star)))
      case('UHZ')
         freq_scale = 1.E6_WP/(TWOPI*SQRT(bc%R_star**3/(bc%G*bc%M_star)))
      case('ACOUSTIC_CUTOFF')
         call eval_cutoffs(bc, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
         freq_scale = 1._WP/omega_cutoff_hi
      case('GRAVITY_CUTOFF')
         call eval_cutoffs(bc, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
         freq_scale = 1._WP/omega_cutoff_lo
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function evol_freq_scale

!****

    function poly_freq_scale (freq_units) result (freq_scale)

      character(LEN=*), intent(in) :: freq_units
      real(WP)                     :: freq_scale

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case (freq_units)
      case ('NONE')
         freq_scale = 1._WP
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function poly_freq_scale

!****

    function hom_freq_scale (freq_units) result (freq_scale)

      character(LEN=*), intent(in) :: freq_units
      real(WP)                     :: freq_scale

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case (freq_units)
      case ('NONE')
         freq_scale = 1._WP
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function hom_freq_scale

  end function freq_scale

end module gyre_units
