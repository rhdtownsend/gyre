! Module   : gyre_freq
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
  use gyre_hom_model
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_poly_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface omega_inert
     module procedure omega_inert_r_
     module procedure omega_inert_c_
  end interface omega_inert

  interface omega_corot
     module procedure omega_corot_r_
     module procedure omega_corot_c_
  end interface omega_corot

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

  public :: omega_inert
  public :: omega_corot
  public :: omega_from_freq
  public :: freq_from_omega
  public :: eval_cutoff_freqs

  ! Procedures

contains

  $define $OMEGA_INERT $sub

  $local $T $1
  $local $TYPE $2

  function omega_inert_${T}_ (omega_c, Omega_rot, m) result (omega)

    $TYPE(WP), intent(in) :: omega_c
    real(WP), intent(in)  :: Omega_rot
    integer, intent(in)   :: m
    $TYPE(WP)             :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c + m*Omega_rot

    ! Finish

    return

  end function omega_inert_${T}_

  $endsub

  $OMEGA_INERT(r,real)
  $OMEGA_INERT(c,complex)

  !****

  $define $OMEGA_COROT $sub

  $local $T $1
  $local $TYPE $2

  function omega_corot_${T}_ (omega, Omega_rot, m) result (omega_c)

    $TYPE(WP), intent(in) :: omega
    real(WP), intent(in)  :: Omega_rot
    integer, intent(in)   :: m
    $TYPE(WP)             :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = omega - m*Omega_rot

    ! Finish

    return

  end function omega_corot_${T}_

  $endsub

  $OMEGA_COROT(r,real)
  $OMEGA_COROT(c,complex)

  !****

  $define $OMEGA_FROM_FREQ $sub

  $local $T $1
  $local $TYPE $2

  function omega_from_freq_${T}_ (freq, ml, pt_i, pt_o, freq_units, freq_frame, md_p, os_p) result (omega)

    $TYPE(WP), intent(in)               :: freq
    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    $TYPE(WP)                           :: omega

    real(WP)  :: Omega_rot
    $TYPE(WP) :: omega_l
    real(WP)  :: omega_cutoff_lo
    real(WP)  :: omega_cutoff_hi

    ! Calculate the dimensionless inertial-frame frequency omega from
    ! the dimensioned local-frame frequency freq

    ! Calculate the dimensionless frequency in the local frame

    select type (ml)

    class is (evol_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case ('HZ')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case ('UHZ')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))/1E6_WP
       case ('RAD_PER_SEC')
          omega_l = freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case ('CYC_PER_DAY')
          omega_l = TWOPI*freq*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))/86400._WP
       case ('ACOUSTIC_DELTA')
          omega_l = TWOPI*freq*ml%Delta_p(pt_i%x, pt_o%x)
       case ('GRAVITY_DELTA')
          omega_l = TWOPI*freq*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))
       case ('UPPER_DELTA')
          omega_l = TWOPI*freq*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('LOWER_DELTA')
          omega_l = TWOPI*freq*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, pt_o, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_hi
       case ('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, pt_o, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_lo
       case ('ROSSBY_I')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1))
       case ('ROSSBY_O')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1))
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case ('ACOUSTIC_DELTA')
          omega_l = TWOPI*freq*ml%Delta_p(pt_i%x, pt_o%x)
       case ('GRAVITY_DELTA')
          omega_l = TWOPI*freq*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))
       case ('UPPER_DELTA')
          omega_l = TWOPI*freq*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('LOWER_DELTA')
          omega_l = TWOPI*freq*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('ROSSBY_I')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1))
       case ('ROSSBY_O')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1))
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (hom_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case ('ACOUSTIC_DELTA')
          omega_l = TWOPI*freq*ml%Delta_p(pt_i%x, pt_o%x)
       case ('GRAVITY_DELTA')
          omega_l = TWOPI*freq*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))
       case ('UPPER_DELTA')
          omega_l = TWOPI*freq*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('LOWER_DELTA')
          omega_l = TWOPI*freq*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('ROSSBY_I')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1))
       case ('ROSSBY_O')
          omega_l = -freq*2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1))
       case default
          $ABORT(Invalid freq_units)
       end select

    class default

       $ABORT(Invalid ml type)

    end select

    ! Now convert to the inertial frame

    select case (freq_frame)
    case ('INERTIAL')
       omega = omega_l
    case ('COROT_I')
       Omega_rot = ml%coeff(I_OMEGA_ROT, pt_i)
       omega = omega_inert(omega_l, Omega_rot, md_p%m)
    case ('COROT_O')
       Omega_rot = ml%coeff(I_OMEGA_ROT, pt_o)
       omega = omega_inert(omega_l, Omega_rot, md_p%m)
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

  function freq_from_omega_${T}_ (omega, ml, pt_i, pt_o, freq_units, freq_frame, md_p, os_p) result (freq)

    $TYPE(WP), intent(in)               :: omega
    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    $TYPE(WP)                           :: freq

    real(WP)  :: Omega_rot
    $TYPE(WP) :: omega_l
    real(WP)  :: omega_cutoff_lo
    real(WP)  :: omega_cutoff_hi

    ! Calculate the dimensioned local-frame frequency freq from the
    ! dimensionless inertial-frame frequency omega

    ! Convert from the inertial frame

    select case (freq_frame)
    case ('INERTIAL')
       omega_l = omega
    case ('COROT_I')
       Omega_rot = ml%coeff(I_OMEGA_ROT, pt_i)
       omega_l = omega_corot(omega, Omega_rot, md_p%m)
    case ('COROT_O')
       Omega_rot = ml%coeff(I_OMEGA_ROT, pt_o)
       omega_l = omega_corot(omega, Omega_rot, md_p%m)
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
       case ('RAD_PER_SEC')
          freq = omega_l/SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       case ('CYC_PER_DAY')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*86400._WP
       case ('ACOUSTIC_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_p(pt_i%x, pt_o%x))
       case ('GRAVITY_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('UPPER_DELTA')
          freq = omega_l/(TWOPI*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('LOWER_DELTA')
          freq = omega_l/(TWOPI*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, pt_o, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, pt_o, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_lo
       case ('ROSSBY_I')
          freq = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1)))
       case ('ROSSBY_O')
          omega_l = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1)))
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          freq = omega_l
       case ('ACOUSTIC_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_p(pt_i%x, pt_o%x))
       case ('GRAVITY_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('UPPER_DELTA')
          freq = omega_l/(TWOPI*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('LOWER_DELTA')
          freq = omega_l/(TWOPI*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('ROSSBY_I')
          freq = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1)))
       case ('ROSSBY_O')
          omega_l = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1)))
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (hom_model_t)

       select case (freq_units)
       case ('NONE')
          freq = omega_l
       case ('ACOUSTIC_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_p(pt_i%x, pt_o%x))
       case ('GRAVITY_DELTA')
          freq = omega_l/(TWOPI*ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP)))
       case ('UPPER_DELTA')
          freq = omega_l/(TWOPI*MAX(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('LOWER_DELTA')
          freq = omega_l/(TWOPI*MIN(ml%Delta_p(pt_i%x, pt_o%x), ml%Delta_g(pt_i%x, pt_o%x, md_p%l*(md_p%l+1._WP))))
       case ('ROSSBY_I')
          freq = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_i)/(md_p%l*(md_p%l+1)))
       case ('ROSSBY_O')
          omega_l = -omega_l/(2._WP*md_p%m*ml%coeff(I_OMEGA_ROT, pt_o)/(md_p%l*(md_p%l+1)))
       case default
          $ABORT(Invalid freq_units)
       end select

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

  subroutine eval_cutoff_freqs (ml, pt, md_p, os_p, omega_cutoff_lo, omega_cutoff_hi)

    class(model_t), intent(in)   :: ml
    type(point_t), intent(in)    :: pt
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    real(WP), intent(out)        :: omega_cutoff_lo
    real(WP), intent(out)        :: omega_cutoff_hi

    real(WP)      :: V_g
    real(WP)      :: As
    real(WP)      :: U
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

       call eval_atmos_coeffs_unno(ml, pt, V_g, As, U, c_1)
       call eval_atmos_cutoff_freqs(V_g, As, U, c_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case('JCD')

       call eval_atmos_coeffs_jcd(ml, pt, V_g, As, U, c_1)
       call eval_atmos_cutoff_freqs(V_g, As, U, c_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

    case('LUAN')

       call eval_atmos_coeffs_luan(ml, pt, V_g, As, U, c_1)
       call eval_atmos_cutoff_freqs(V_g, As, U, c_1, md_p%l*(md_p%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

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
