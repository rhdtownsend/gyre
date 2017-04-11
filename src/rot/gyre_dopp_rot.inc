! Incfile  : gyre_dopp_rot
! Purpose  : rotational effects (Doppler shift, template)
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_${T}_dopp_rot

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_mode_par
  use gyre_point
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_dopp_rot_t
     private
     class(model_t), pointer :: ml => null()
     real(WP), allocatable   :: Omega_rot(:)
     integer                 :: l
     integer                 :: m
   contains
     private
     procedure, public :: stencil
     procedure, public :: omega
     procedure, public :: omega_c
     procedure, public :: l_i
     procedure, public :: l_e
     procedure, public :: lambda
  end type ${T}_dopp_rot_t

  ! Interfaces

  interface ${T}_dopp_rot_t
     module procedure ${T}_dopp_rot_t_
  end interface ${T}_dopp_rot_t

  ! Access specifiers

  private

  public :: ${T}_dopp_rot_t

  ! Procedures

contains

  function ${T}_dopp_rot_t_ (ml, md_p) result (rt)

    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(${T}_dopp_rot_t)               :: rt

    ! Construct the dopp_rot_t

    rt%ml => ml

    rt%l = md_p%l
    rt%m = md_p%m

    ! Finish

    return

  end function ${T}_dopp_rot_t_

  !****

  subroutine stencil (this, pt)

    class(${T}_dopp_rot_t), intent(inout) :: this
    type(point_t), intent(in)             :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    n_s = SIZE(pt)

    if (ALLOCATED(this%Omega_rot)) deallocate(this%Omega_rot)
    allocate(this%Omega_rot(n_s))

    do i = 1, n_s
       this%Omega_rot(i) = this%ml%coeff(I_OMEGA_ROT, pt(i))
    end do

    ! Finish

    return

  end subroutine stencil

  !****

  function omega (this, i, omega_c)

    class(${T}_dopp_rot_t), intent(in) :: this
    integer, intent(in)                :: i
    $TYPE(WP), intent(in)              :: omega_c
    $TYPE(WP)                          :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c + this%m*this%Omega_rot(i)

    ! Finish

    return

  end function omega

  !****

  function omega_c (this, i, omega)

    class(${T}_dopp_rot_t), intent(in) :: this
    integer, intent(in)                :: i
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = omega - this%m*this%Omega_rot(i)

    ! Finish

    return

  end function omega_c

  !****

  function l_i (this, omega)

    class(${T}_dopp_rot_t), intent(in) :: this
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_i

    ! Evaluate the effective harmonic degree at the model inner
    ! boundary

    l_i = this%l

    ! Finish

    return

  end function l_i

  !****

  function l_e (this, i, omega)

    class(${T}_dopp_rot_t), intent(in) :: this
    integer, intent(in)                :: i
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_e

    ! Evaluate the effective harmonic degree

    l_e = this%l

    ! Finish

    return
    
  end function l_e

  !****

  function lambda (this, i, omega)

    class(${T}_dopp_rot_t), intent(in) :: this
    integer, intent(in)                :: i
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: lambda

    ! Evaluate the angular eigenvalue

    lambda = this%l*(this%l+1)

    ! Finish

    return

  end function lambda

end module gyre_${T}_dopp_rot

$endsub
