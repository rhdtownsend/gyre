! Module   : gyre_freq_context
! Purpose  : frequency transformation routines (context-dependent)
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

module gyre_freq_context

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_constants
  use gyre_context
  use gyre_evol_model
  use gyre_freq_model
  use gyre_math
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface freq_scale
     module procedure freq_scale_cx_
  end interface freq_scale

  interface freq_shift
     module procedure freq_shift_cx_
  end interface freq_shift

  ! Access specifiers

  private

  public :: freq_scale
  public :: freq_shift

  ! Procedures

contains

  function freq_scale_cx_ (units, cx, md_p, os_p) result (scale)

    character(*), intent(in)     :: units
    class(context_t), intent(in) :: cx
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    real(WP)                     :: scale

    class(model_t), pointer :: ml => null()
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o

    ! Evaluate the multiplicative scale for converting dimensioned
    ! frequencies (as specified by units) to dimensionless angular
    ! frequencies

    ml => cx%model()
    
    pt_i = cx%point_i()
    pt_o = cx%point_o()

    select case (units)
    case ('ACOUSTIC_DELTA')
       scale = TWOPI*ml%Delta_p(pt_i%x, pt_o%x)
    case ('GRAVITY_DELTA')
       scale = TWOPI*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))
    case ('UPPER_DELTA')
       scale = TWOPI*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
    case ('LOWER_DELTA')
       scale = TWOPI*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
    case ('ROSSBY_I')
       scale = -2._WP*md_p%m*cx%Omega_rot(pt_i)/(md_p%l*(md_p%l+1))
    case ('ROSSBY_O')
       scale = -2._WP*md_p%m*cx%Omega_rot(pt_o)/(md_p%l*(md_p%l+1))
    case default
       select type (ml)
       class is (evol_model_t)
          scale = freq_scale_cx_evol_(units, ml, pt_o, md_p, os_p)
       class default
          scale = freq_scale(units, ml)
       end select
    end select

    ! Finish

    return

  end function freq_scale_cx_

  !****

  function freq_scale_cx_evol_ (units, ml, pt, md_p, os_p) result (scale)

    character(*), intent(in)       :: units
    type(evol_model_t), intent(in) :: ml
    type(point_t), intent(in)      :: pt
    type(mode_par_t), intent(in)   :: md_p
    type(osc_par_t), intent(in)    :: os_p
    real(WP)                       :: scale

    real(WP) :: omega_cutoff_lo
    real(WP) :: omega_cutoff_hi

    select case (units)
    case ('ACOUSTIC_CUTOFF')
       call eval_cutoff_freqs(ml, pt, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
       scale = omega_cutoff_hi
    case ('GRAVITY_CUTOFF')
       call eval_cutoff_freqs(ml, pt, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
       scale = omega_cutoff_lo
    case default
       scale = freq_scale(units, ml)
    end select

    ! Finish

    return

  end function freq_scale_cx_evol_

  !****

  function freq_shift_cx_ (frame, cx, md_p) result (shift)

    character(*), intent(in)     :: frame
    class(context_t), intent(in) :: cx
    type(mode_par_t), intent(in) :: md_p
    real(WP)                     :: shift
    
    class(model_t), pointer :: ml => null()
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o

    ! Evaluate the additive shift for converting dimensionless angular
    ! frequencies from a co-rotating frame (as specified by frame) to
    ! the inertial frame
    
    ml => cx%model()
    
    pt_i = cx%point_i()
    pt_o = cx%point_o()

    ! Now convert to the inertial frame

    select case (frame)
    case ('INERTIAL')
       shift = 0._WP
    case ('COROT_I')
       shift = md_p%m*cx%Omega_rot(pt_i)
    case ('COROT_O')
       shift = md_p%m*cx%Omega_rot(pt_o)
    case default
       $ABORT(Invalid frequency frame)
    end select

    ! Finish

    return

  end function freq_shift_cx_
    
  !****

  subroutine eval_cutoff_freqs (ml, pt, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)

    class(model_t), intent(in)   :: ml
    type(point_t), intent(in)    :: pt
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    real(WP), intent(out)        :: omega_cutoff_lo
    real(WP), intent(out)        :: omega_cutoff_hi

    character(LEN(os_p%outer_bound)) :: outer_bound
    real(WP)                         :: V
    real(WP)                         :: As
    real(WP)                         :: c_1
    real(WP)                         :: Gamma_1
    logical, save                    :: warned = .FALSE.

    ! Evaluate the cutoff frequencies

    if (os_p%outer_bound_for_cutoff /= '') then
       outer_bound = os_p%outer_bound_for_cutoff
    else
       outer_bound = os_p%outer_bound
    endif

    select case (os_p%outer_bound)

    case ('VACUUM')

       $ABORT(Cutoff frequencies cannot be evaluated for VACUUM outer boundary condition)

    case ('DZIEM')

       $ABORT(Cutoff frequencies cannot be evaluated for DZIEM outer boundary condition)

    case ('ISOTHERMAL')

       call eval_atmos_coeffs_isothrm(ml, pt, V, As, c_1, Gamma_1)
       call eval_atmos_cutoff_freqs(V, As, c_1, Gamma_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case ('UNNO')

       call eval_atmos_coeffs_unno(ml, pt, V, As, c_1, Gamma_1)
       call eval_atmos_cutoff_freqs(V, As, c_1, Gamma_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case ('JCD')

       call eval_atmos_coeffs_isothrm(ml, pt, V, As, c_1, Gamma_1)
       call eval_atmos_cutoff_freqs(V, As, c_1, Gamma_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case default

       $ABORT(Invalid outer_bound)

    end select

    if (.not. warned) then
       $WARN(WARNING: Cutoff frequencies do not account for rotation effects)
       warned = .TRUE.
    endif

    ! Finish

    return

  end subroutine eval_cutoff_freqs

end module gyre_freq_context
