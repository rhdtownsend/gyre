! Incfile  : gyre_null_rot
! Purpose  : rotational effects (null, template)
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

module gyre_${T}_null_rot

  ! Uses

  use core_kinds

  use gyre_modepar
  use gyre_rot
 
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_null_rot_t
     private
   contains
     private
     procedure, public :: omega => omega_
     procedure         :: omega_c_1_
     procedure         :: omega_c_v_
     procedure         :: lambda_1_
     procedure         :: lambda_v_
  end type ${T}_null_rot_t

  ! Interfaces

  interface ${T}_null_rot_t
     module procedure ${T}_null_rot_t_
  end interface ${T}_null_rot_t

  ! Access specifiers

  private

  public :: ${T}_null_rot_t

  ! Procedures

contains

  function ${T}_null_rot_t_ (mp) result (rt)

    type(modepar_t), intent(in) :: mp
    type(${T}_null_rot_t)       :: rt

    ! Construct the null_rot_t

    rt%mp = mp

    ! Finish

    return

  end function ${T}_null_rot_t_

!****

  function omega_ (this, x, omega_c) result (omega)

    class(${T}_null_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega_c
    $TYPE(WP)                          :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c

    ! Finish

    return

  end function omega_

!****

  function omega_c_1_ (this, x, omega) result (omega_c)

    class(${T}_null_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency (scalar)

    omega_c = omega

    ! Finish

    return

  end function omega_c_1_

!****

  function omega_c_v_ (this, x, omega) result (omega_c)

    class(${T}_null_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x(:)
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c(SIZE(x))

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency (vector)

    omega_c = omega

    ! Finish

    return

  end function omega_c_v_

!****

  function lambda_1_ (this, x, omega) result (lambda)

    class(${T}_null_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: lambda

    ! Evaluate the angular eigenvalue (scalar)

    lambda = this%mp%l*(this%mp%l+1)

    ! Finish

    return

  end function lambda_1_

!****

  function lambda_v_ (this, x, omega) result (lambda)

    class(${T}_null_rot_t), intent(in) :: this
    real(WP), intent(in)               :: x(:)
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: lambda(SIZE(x))

    ! Evaluate the angular eigenvalue (vector)

    lambda = this%mp%l*(this%mp%l+1)

    ! Finish

    return

  end function lambda_v_

end module gyre_${T}_null_rot

$endsub
