! Incfile  : gyre_rot
! Purpose  : rotational effects (template)
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

module gyre_${T}_rot

  ! Uses

  use core_kinds

  use gyre_modepar

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_rot_t
     private
     type(modepar_t), public :: mp
   contains
     private
     procedure(omega_), deferred, public :: omega
     procedure(omega_c_1_), deferred     :: omega_c_1_
     procedure(omega_c_v_), deferred     :: omega_c_v_
     generic, public                     :: omega_c => omega_c_1_, omega_c_v_
     procedure                           :: l_e_1_
     procedure                           :: l_e_v_
     generic, public                     :: l_e => l_e_1_, l_e_v_
     procedure, public                   :: l_0 => l_0_
     procedure(lambda_1_), deferred      :: lambda_1_
     procedure(lambda_v_), deferred      :: lambda_v_
     generic, public                     :: lambda => lambda_1_, lambda_v_
  end type ${T}_rot_t

  ! Interfaces

  abstract interface

     function omega_ (this, x, omega_c) result (omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega_c
       $TYPE(WP)                     :: omega
     end function omega_

     function omega_c_1_ (this, x, omega) result (omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c
     end function omega_c_1_

     function omega_c_v_ (this, x, omega) result (omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c(SIZE(x))
     end function omega_c_v_

     function lambda_1_ (this, x, omega) result (lambda)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda
     end function lambda_1_
     
     function lambda_v_ (this, x, omega) result (lambda)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda(SIZE(x))
     end function lambda_v_
     
  end interface

  ! Access specifiers

  private

  public :: ${T}_rot_t

contains

  function l_e_1_ (this, x, omega) result (l_e)

    class(${T}_rot_t), intent(in) :: this
    real(WP), intent(in)          :: x
    $TYPE(WP), intent(in)         :: omega
    $TYPE(WP)                     :: l_e

    ! Evaluate the effective harmonic degree (scalar)

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda(x, omega)))

    ! Finish

    return

  end function l_e_1_

!****

  function l_e_v_ (this, x, omega) result (l_e)

    class(${T}_rot_t), intent(in) :: this
    real(WP), intent(in)          :: x(:)
    $TYPE(WP), intent(in)         :: omega
    $TYPE(WP)                     :: l_e(SIZE(x))

    ! Evaluate the effective harmonic degree (vector)

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda(x, omega)))

    ! Finish

    return

  end function l_e_v_

!****

  function l_0_ (this, omega) result (l_0)

    class(${T}_rot_t), intent(in) :: this
    $TYPE(WP), intent(in)         :: omega
    $TYPE(WP)                     :: l_0

    ! Evaluate the effective harmonic degree at the origin

    l_0 = this%l_e(0._WP, omega)

    ! Finish

    return

  end function l_0_

end module gyre_${T}_rot

$endsub
