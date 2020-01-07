! Module   : gyre_freq_frame
! Purpose  : frequency frame transformation routines
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

module gyre_freq_frame

  ! Uses

  use core_kinds

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

  ! Access specifiers

  private

  public :: omega_inert
  public :: omega_corot

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

end module gyre_freq_frame
