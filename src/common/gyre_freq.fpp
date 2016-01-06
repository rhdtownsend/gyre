! Program  : gyre_freq
! Purpose  : frequency transformation routines
!
! Copyright 2015 Rich Townsend
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

module gyre_freq

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_constants
  use gyre_evol_model
!  use gyre_hom_model
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_poly_model
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface omega_from_freq
     module procedure omega_from_freq_r_
     module procedure omega_from_freq_c_
  end interface omega_from_freq

  interface freq_from_omega
     module procedure freq_from_omega_r_
     module procedure freq_from_omega_c_
  end interface freq_from_omega

  ! Access specifiers

  private

  public :: omega_from_freq
  public :: freq_from_omega
  public :: eval_cutoff_freqs

  ! Procedures

contains

  $define $OMEGA_FROM_FREQ $sub

  $local $T $1
  $local $TYPE $2

  function omega_from_freq_${T}_ (freq, ml, freq_units, freq_frame, md_p, os_p) result (omega)

    $TYPE(WP), intent(in)               :: freq
    class(model_t), pointer, intent(in) :: ml
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    $TYPE(WP)                           :: omega

    $TYPE(WP)                      :: omega_l
    class(${T}_rot_t), allocatable :: rt
    real(WP)                       :: omega_cutoff_lo
    real(WP)                       :: omega_cutoff_hi

    ! Calculate the dimensionless inertial-frame frequency omega from
    ! the dimensioned local-frame frequency freq

    ! First calculate the dimensionless frequency in the local frame

    select type (ml)

    class is (evol_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case ('HZ')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case ('UHZ')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))/1E6_WP
       case ('PER_DAY')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))/86400._WP
       case ('ACOUSTIC_DELTA')
          omega_l = TWOPI*freq*ml%delta_p()*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case ('GRAVITY_DELTA')
          omega_l = TWOPI*freq*ml%delta_g()*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case default
          $ABORT(Invalid freq_units)
       end select

   ! class is (hom_model_t)

   !     select case (freq_units)
   !     case ('NONE')
   !        omega_l = freq
   !     case default
   !       $ABORT(Invalid freq_units)
   !    end select

    class default

       $ABORT(Invalid ml type)

    end select

    ! Now convert to the inertial frame

    allocate(rt, SOURCE=${T}_rot_t(ml, md_p, os_p))

    select case (freq_frame)
    case ('INERTIAL')
       omega = omega_l
    case ('COROT_I')
       omega = rt%omega(1, ml%x_i(1), omega_l)
    case ('COROT_O')
       omega = rt%omega(ml%n_s, ml%x_o(ml%n_s), omega_l)
    case default
       $ABORT(Invalid freq_frame)
    end select

    ! Finish

    return

  end function omega_from_freq_${T}_

  $endsub

  $OMEGA_FROM_FREQ(r,real)
  $OMEGA_FROM_FREQ(c,complex)

!****

  $define $FREQ_FROM_OMEGA $sub

  $local $T $1
  $local $TYPE $2

  function freq_from_omega_${T}_ (omega, ml, freq_units, freq_frame, md_p, os_p) result (freq)

    $TYPE(WP), intent(in)               :: omega
    class(model_t), pointer, intent(in) :: ml
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    $TYPE(WP)                           :: freq

    $TYPE(WP)                      :: omega_l
    class(${T}_rot_t), allocatable :: rt
    real(WP)                       :: omega_cutoff_lo
    real(WP)                       :: omega_cutoff_hi

    ! Calculate the dimensioned local-frame frequency freq from the
    ! dimensionless inertial-frame frequency omega

    ! First convert from the inertial frame

    allocate(rt, SOURCE=${T}_rot_t(ml, md_p, os_p))

    select case (freq_frame)
    case ('INERTIAL')
       omega_l = omega
    case ('COROT_I')
       omega_l = rt%omega_c(1, ml%x_i(1), omega)
    case ('COROT_O')
       omega_l = rt%omega_c(ml%n_s, ml%x_o(ml%n_s), omega)
    case default
       $ABORT(Invalid freq_frame)
    end select

    ! Now calculate the dimensioned frequency in the local frame

    select type (ml)

    class is (evol_model_t)

       select case(freq_units)
       case ('NONE')
          freq = omega_l
       case ('HZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case ('UHZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*1E6_WP
       case ('PER_DAY')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*86400._WP
       case ('ACOUSTIC_DELTA')
          freq = omega_l/(TWOPI*ml%delta_p()*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case ('GRAVITY_DELTA')
          freq = omega_l/(TWOPI*ml%delta_g()*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case ('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          freq = omega_l
       case default
          $ABORT(Invalid freq_units)
       end select

    ! class is (hom_model_t)

    !    select case (freq_units)
    !    case ('NONE')
    !       freq = omega_l
    !    case default
    !      $ABORT(Invalid freq_units)
    !   end select

    class default

       $ABORT(Invalid ml type)

    end select

    ! Finish

    return

  end function freq_from_omega_${T}_

  $endsub

  $FREQ_FROM_OMEGA(r,real)
  $FREQ_FROM_OMEGA(c,complex)

!****

  subroutine eval_cutoff_freqs (ml, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)

    class(model_t), intent(in)   :: ml
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    real(WP), intent(out)        :: omega_cutoff_lo
    real(WP), intent(out)        :: omega_cutoff_hi
    
    real(WP)      :: V_g
    real(WP)      :: As
    real(WP)      :: c_1
    logical, save :: warned = .FALSE.

    ! Evaluate the cutoff frequencies

    select case (os_p%outer_bound)

    case ('ZERO')

       omega_cutoff_lo = 0._WP
       omega_cutoff_hi = HUGE(0._WP)

    case ('DZIEM')

       omega_cutoff_lo = 0._WP
       omega_cutoff_hi = HUGE(0._WP)

    case ('UNNO')

       call eval_atmos_coeffs_unno(ml, V_g, As, c_1)
       call eval_atmos_cutoff_freqs(V_g, As, c_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case('JCD')

       call eval_atmos_coeffs_jcd(ml, V_g, As, c_1)
       call eval_atmos_cutoff_freqs(V_g, As, c_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

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

end module gyre_freq
