! Incfile  : gyre_trad_rot
! Purpose  : rotational effects (traditional approximation, template)
!
! Copyright 2013-2014 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_trad_rot

  ! Uses

  use core_kinds

  use gyre_modepar
  use gyre_model
  use gyre_trad
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_trad_rot_t
     private
     class(model_t), pointer :: ml => null()
     type(modepar_t)         :: mp
   contains
     private
     procedure, public :: l_e => l_e_
     procedure, public :: l_0 => l_0_
     procedure, public :: omega => omega_
     procedure, public :: omega_c => omega_c_
  end type ${T}_trad_rot_t

  ! Interfaces

  interface ${T}_trad_rot_t
     module procedure ${T}_trad_rot_t_
  end interface ${T}_trad_rot_t

  ! Access specifiers

  private

  public :: ${T}_trad_rot_t

  ! Procedures

contains

  function ${T}_trad_rot_t_ (ml, mp) result (wv)

    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(${T}_trad_rot_t)               :: wv

    ! Construct the trad_rot_t

    wv%ml => ml
    wv%mp = mp

    ! Finish

    return

  end function ${T}_trad_rot_t_

!****

  function l_e_ (this, x, omega) result (l_e)

    class(${T}_trad_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_e

    $TYPE(WP) :: nu
    $TYPE(WP) :: lambda

    ! Evaluate the effective harmonic degree

    nu = 2._WP*this%ml%Omega_rot(x)/this%omega_c(x, omega)

    lambda = trad_lambda(nu, this%mp%l, this%mp%m)

    Lx_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4.*lambda))

    ! Finish

    return

  end function l_e_

!****

  function l_0_ (this, omega) result (l_0)

    class(${T}_trad_rot_t), intent(in) :: this
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_0

    ! Evaluate the harmonic degree used to scale solutions

    l_0 = this%l_e(0._WP, omega)

    ! Finish

    return

  end function l_0_

!****

  function omega_ (this, x, omega_c) result (omega)

    class(${T}_trad_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega_c
    $TYPE(WP)                          :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c + this%mp%m*this%ml%Omega_rot(x)

    ! Finish

    return

  end function omega_

!****

  function omega_c_ (this, x, omega) result (omega_c)

    class(${T}_trad_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = omega - this%mp%m*this%ml%Omega_rot(x)

    ! Finish

    return

  end function omega_c_

end module gyre_${T}_trad_rot

$endsub
