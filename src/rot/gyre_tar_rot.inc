! Incfile  : gyre_tar_rot
! Purpose  : rotational effects (traditional approximation of rotation, template)
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

module gyre_${T}_tar_rot

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_mode_par
  use gyre_point
  use gyre_tar
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_tar_rot_t
     private
     class(model_t), pointer :: ml => null()
     real(WP), allocatable   :: Omega_rot(:)
     real(WP)                :: Omega_rot_i
     integer                 :: l
     integer                 :: m
     logical                 :: rossby
   contains
     private
     procedure, public :: stencil
     procedure, public :: omega
     procedure, public :: omega_c
     procedure         :: omega_c_O_r_
     procedure, public :: l_i
     procedure, public :: l_e
     procedure         :: l_e_O_r_
     procedure, public :: lambda
     procedure         :: lambda_O_r_
  end type ${T}_tar_rot_t

  ! Interfaces

  interface ${T}_tar_rot_t
     module procedure ${T}_tar_rot_t_
  end interface ${T}_tar_rot_t

  ! Access specifiers

  private

  public :: ${T}_tar_rot_t

  ! Procedures

contains

  function ${T}_tar_rot_t_ (ml, pt_i, md_p) result (rt)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(mode_par_t), intent(in)        :: md_p
    type(${T}_tar_rot_t)                :: rt

    ! Construct the tar_rot_t

    rt%ml => ml

    rt%Omega_rot_i = ml%coeff(I_OMEGA_ROT, pt_i)

    rt%l = md_p%l
    rt%m = md_p%m

    rt%rossby = md_p%rossby

    ! Finish

    return

  end function ${T}_tar_rot_t_

  !****

  subroutine stencil (this, pt)

    class(${T}_tar_rot_t), intent(inout) :: this
    type(point_t), intent(in)            :: pt(:)

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

    class(${T}_tar_rot_t), intent(in) :: this
    integer, intent(in)               :: i
    $TYPE(WP), intent(in)             :: omega_c
    $TYPE(WP)                         :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c + this%m*this%Omega_rot(i)

    ! Finish

    return

  end function omega

  !****

  function omega_c (this, i,  omega)

    class(${T}_tar_rot_t), intent(in) :: this
    integer, intent(in)               :: i
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = this%omega_c_O_r_(this%Omega_rot(i), omega)

    ! Finish

    return

  end function omega_c

  !****

  function omega_c_O_r_ (this, Omega_rot, omega) result (omega_c)

    class(${T}_tar_rot_t), intent(in) :: this
    real(WP), intent(in)              :: Omega_rot
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency for the given rotation rate

    omega_c = omega - this%m*Omega_rot

    ! Finish

    return

  end function omega_c_O_r_

  !****

  function l_i (this, omega)

    class(${T}_tar_rot_t), intent(in) :: this
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: l_i

    ! Evaluate the effective harmonic degree at the inner boundary

    l_i = this%l_e_O_r_(this%Omega_rot_i, omega)

    ! Finish

    return

  end function l_i

  !****

  function l_e (this, i, omega)

    class(${T}_tar_rot_t), intent(in) :: this
    integer, intent(in)               :: i
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: l_e

    l_e = this%l_e_O_r_(this%Omega_rot(i), omega)

    ! Finish

    return

  end function l_e

  !****

  function l_e_O_r_ (this, Omega_rot, omega) result (l_e)

    class(${T}_tar_rot_t), intent(in) :: this
    real(WP), intent(in)              :: Omega_rot
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: l_e

    ! Evaluate the effective harmonic degree for the given rotation
    ! rate

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda_O_r_(Omega_rot, omega)))

    ! Finish

    return

  end function l_e_O_r_

  !****

  function lambda (this, i, omega)

    class(${T}_tar_rot_t), intent(in) :: this
    integer, intent(in)               :: i
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: lambda

    ! Evaluate the angular eigenvalue

    lambda = this%lambda_O_r_(this%Omega_rot(i), omega)

    ! Finish

    return

  end function lambda

  !****
  
  function lambda_O_r_ (this, Omega_rot, omega) result (lambda)

    class(${T}_tar_rot_t), intent(in) :: this
    real(WP), intent(in)              :: Omega_rot
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: lambda

    $TYPE(WP) :: nu

    ! Evaluate the angular eigenvalue for the given rotation rate

    nu = 2._WP*Omega_rot/this%omega_c_O_r_(Omega_rot, omega)

    lambda = tar_lambda(nu, this%l, this%m, this%rossby)

    ! Finish

    return

  end function lambda_O_r_

end module gyre_${T}_tar_rot

$endsub
