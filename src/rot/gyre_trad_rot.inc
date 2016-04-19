! Incfile  : gyre_trad_rot
! Purpose  : rotational effects (traditional approximation, template)
!
! Copyright 2013-2015 Rich Townsend
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

  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_point
  use gyre_trad
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_trad_rot_t
     private
     class(model_t), pointer :: ml => null()
     type(point_t)           :: pt_i
     integer                 :: l
     integer                 :: m
   contains
     private
     procedure, public :: omega
     procedure         :: omega_c_1_
     procedure         :: omega_c_v_
     procedure, public :: l_i
     procedure         :: l_e_1_
     procedure         :: l_e_v_
     procedure         :: lambda_1_
     procedure         :: lambda_v_
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

  function ${T}_trad_rot_t_ (ml, gr, md_p) result (rt)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(${T}_trad_rot_t)               :: rt

    ! Construct the trad_rot_t

    rt%ml => ml

    rt%pt_i = gr%pt(1)

    rt%l = md_p%l
    rt%m = md_p%m

    ! Finish

    return

  end function ${T}_trad_rot_t_

  !****

  function omega (this, pt, omega_c)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    $TYPE(WP), intent(in)              :: omega_c
    $TYPE(WP)                          :: omega

    ! Evaluate the inertial-frame frequency from the corotating-frame
    ! frequency

    omega = omega_c + this%m*this%ml%Omega_rot(pt)

    ! Finish

    return

  end function omega

  !****

  function omega_c_1_ (this, pt,  omega) result (omega_c)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = omega - this%m*this%ml%Omega_rot(pt)

    ! Finish

    return

  end function omega_c_1_

  !****

  function omega_c_v_ (this, pt, omega) result (omega_c)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt(:)
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: omega_c(SIZE(pt))

    ! Evaluate the corotating-frame frequency from the inertial-frame
    ! frequency

    omega_c = omega - this%m*this%ml%Omega_rot(pt)

    ! Finish

    return

  end function omega_c_v_

  !****

  function l_i (this, omega)

    class(${T}_trad_rot_t), intent(in) :: this
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_i

    ! Evaluate the effective harmonic degree at the model inner
    ! boundary

    l_i = this%l_e(this%pt_i, omega)

    ! Finish

    return

  end function l_i

  !****

  function l_e_1_ (this, pt, omega) result (l_e)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_e

    ! Evaluate the effective harmonic degree

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda(pt, omega)))

    ! Finish

    return

  end function l_e_1_

  !****

  function l_e_v_ (this, pt, omega) result (l_e)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt(:)
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: l_e(SIZE(pt))

    ! Evaluate the effective harmonic degree

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda(pt, omega)))

    ! Finish

    return

  end function l_e_v_

  !****

  function lambda_1_ (this, pt, omega) result (lambda)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: lambda

    $TYPE(WP) :: nu

    ! Evaluate the angular eigenvalue

    nu = 2._WP*this%ml%Omega_rot(pt)/this%omega_c(pt, omega)

    lambda = trad_lambda(nu, this%l, this%m)

    ! Finish

    return

  end function lambda_1_

  !****

  function lambda_v_ (this, pt, omega) result (lambda)

    class(${T}_trad_rot_t), intent(in) :: this
    type(point_t), intent(in)          :: pt(:)
    $TYPE(WP), intent(in)              :: omega
    $TYPE(WP)                          :: lambda(SIZE(pt))

    integer :: j
    
    ! Evaluate the angular eigenvalue

    do j = 1, SIZE(pt)
       lambda(j) = this%lambda(pt(j), omega)
    end do

    ! Finish

    return

  end function lambda_v_

end module gyre_${T}_trad_rot

$endsub
