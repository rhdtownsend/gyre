! Module   : gyre_freq
! Purpose  : frequency transformation routines
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

module gyre_freq

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_context
  use gyre_freq_model
  use gyre_freq_context
  use gyre_mode_par
  use gyre_osc_par

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

  ! Procedures

contains

  $define $OMEGA_FROM_FREQ $sub

  $local $T $1
  $local $TYPE $2

  function omega_from_freq_${T}_ (freq, cx, freq_units, freq_frame, md_p, os_p) result (omega)

    $TYPE(WP), intent(in)        :: freq
    class(context_t), intent(in) :: cx
    character(*), intent(in)     :: freq_units
    character(*), intent(in)     :: freq_frame
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    $TYPE(WP)                    :: omega

    omega = freq_scale(freq_units, cx, md_p, os_p)*freq + freq_shift(freq_frame, cx, md_p)

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

  function freq_from_omega_${T}_ (omega, cx, freq_units, freq_frame, md_p, os_p) result (freq)

    $TYPE(WP), intent(in)        :: omega
    class(context_t), intent(in) :: cx
    character(*), intent(in)     :: freq_units
    character(*), intent(in)     :: freq_frame
    type(mode_par_t), intent(in) :: md_p
    type(osc_par_t), intent(in)  :: os_p
    $TYPE(WP)                    :: freq

    ! Calculate the dimensioned local-frame frequency freq from the
    ! dimensionless inertial-frame frequency omega

    freq = (omega - freq_shift(freq_frame, cx, md_p))/freq_scale(freq_units, cx, md_p, os_p)

    ! Finish

    return

  end function freq_from_omega_${T}_

  $endsub

  $FREQ_FROM_OMEGA(r,real)
  $FREQ_FROM_OMEGA(c,complex)

end module gyre_freq
