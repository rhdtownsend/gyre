! Module   : gyre_ad_share
! Purpose  : adiabatic shared data
!
! Copyright 2017 Rich Townsend
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

module gyre_ad_share

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_share_t
     private
     class(model_t), pointer, public :: ml => null()
     class(r_rot_t), allocatable     :: rt
     real(WP)                        :: Omega_rot_i
   contains
     private
     procedure, public :: omega_c
     procedure, public :: lambda
     procedure, public :: l_e
     procedure, public :: l_i
  end type ad_share_t

  ! Interfaces

  interface ad_share_t
     module procedure ad_share_t_
  end interface ad_share_t

  ! Access specifiers

  private

  public :: ad_share_t

  ! Procedures

contains

  function ad_share_t_ (ml, pt_i, md_p, os_p) result (sh)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_share_t)                    :: sh

    ! Construct the ad_share_t

    sh%ml => ml
 
    allocate(sh%rt, SOURCE=r_rot_t(md_p, os_p))

    sh%Omega_rot_i = ml%coeff(I_OMEGA_ROT, pt_i)

    ! Finish

    return

  end function ad_share_t_

  !****

  function omega_c (this, Omega_rot, omega)

    class(ad_share_t), intent(in) :: this
    real(WP), intent(in)          :: Omega_rot
    real(WP), intent(in)          :: omega
    real(WP)                      :: omega_c

    ! Evaluate the corotating-frame frequency

    omega_c = this%rt%omega_c(Omega_rot, omega)

    ! Finish

    return

  end function omega_c

  !****

  function lambda (this, Omega_rot, omega)

    class(ad_share_t), intent(in) :: this
    real(WP), intent(in)          :: Omega_rot
    real(WP), intent(in)          :: omega
    real(WP)                      :: lambda

    ! Evaluate the angular eigenvalue

    lambda = this%rt%lambda(Omega_rot, omega)

    ! Finish

    return

  end function lambda

  !****

  function l_e (this, Omega_rot, omega)

    class(ad_share_t), intent(in) :: this
    real(WP), intent(in)          :: Omega_rot
    real(WP), intent(in)          :: omega
    real(WP)                      :: l_e

    ! Evaluate the effective harmonic degree

    l_e = this%rt%l_e(Omega_rot, omega)

    ! Finish

    return

  end function l_e

  !****

  function l_i (this, omega)

    class(ad_share_t), intent(in) :: this
    real(WP), intent(in)         :: omega
    real(WP)                     :: l_i

    ! Evaluate the effective harmonic degree at the inner boundary

    l_i = this%l_e(this%Omega_rot_i, omega)

    ! Finish

    return

  end function l_i

end module gyre_ad_share
