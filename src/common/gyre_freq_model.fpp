! Module   : gyre_freq_model
! Purpose  : frequency transformation routines (model-dependent)
!
! Copyright 2015-2020 Rich Townsend & The GYRE Team
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

module gyre_freq_model

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_evol_model
  use gyre_math
  use gyre_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface freq_scale
     module procedure freq_scale_ml_
  end interface freq_scale

  ! Access specifiers

  private

  public :: freq_scale

  ! Procedures

contains

  function freq_scale_ml_ (units, ml) result (scale)

    character(*), intent(in)   :: units
    class(model_t), intent(in) :: ml
    real(WP)                   :: scale

    ! Evaluate the multiplicative scale for converting dimensioned
    ! frequencies (as specified by units) to dimensionless angular
    ! frequencies

    select case (units)
    case ('NONE')
       scale = 1._WP
    case default
       select type (ml)
       class is (evol_model_t)
          scale = freq_scale_ml_evol_(units, ml)
       class default
          $ABORT(Invalid frequency units)
       end select
    end select

    ! Finish

    return

  end function freq_scale_ml_

  !****

  function freq_scale_ml_evol_ (units, ml) result (scale)

    character(*), intent(in)       :: units
    type(evol_model_t), intent(in) :: ml
    real(WP)                       :: scale

    select case (units)
    case ('HZ')
       scale = TWOPI*sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))
    case ('UHZ')
       scale = TWOPI*sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))/1E6_WP
    case ('RAD_PER_SEC')
       scale = sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))
    case ('CYC_PER_DAY')
       scale = TWOPI*sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))/86400._WP
    case default
       $ABORT(Invalid frequency units)
    end select

    ! Finish

    return

  end function freq_scale_ml_evol_
  
end module gyre_freq_model
