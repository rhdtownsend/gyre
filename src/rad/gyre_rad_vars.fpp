! Module   : gyre_rad_vars
! Purpose  : adiabatic radial variables transformations
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

module gyre_rad_vars

  ! Uses

  use core_kinds

  use gyre_linalg
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GYRE_SET = 0
  integer, parameter :: DZIEM_SET = 1
  integer, parameter :: JCD_SET = 2
  integer, parameter :: MIX_SET = 3
  integer, parameter :: LAGP_SET = 4

  integer, parameter :: J_V_2 = 1
  integer, parameter :: J_DV_2 = 2
  integer, parameter :: J_C_1 = 3
  integer, parameter :: J_DC_1 = 4

  integer, parameter :: J_LAST = J_DC_1

  ! Derived-type definitions

  type :: rad_vars_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     real(WP), allocatable       :: coeffs(:,:)
     integer                     :: set
   contains
     private
     procedure, public :: stencil
     procedure, public :: G
     procedure         :: G_jcd_
     procedure         :: G_lagp_
     procedure, public :: H
     procedure         :: H_jcd_
     procedure         :: H_lagp_
     procedure, public :: dH
     procedure         :: dH_jcd_
     procedure         :: dH_lagp_
  end type rad_vars_t

  ! Interfaces

  interface rad_vars_t
     module procedure rad_vars_t_
  end interface rad_vars_t

  ! Access specifiers

  private

  public :: rad_vars_t

  ! Procedures

contains

  function rad_vars_t_ (ml, pt_i, md_p, os_p) result (vr)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_vars_t)                    :: vr

    ! Construct the rad_vars_t

    vr%ml => ml

    allocate(vr%rt, SOURCE=r_rot_t(ml, pt_i, md_p, os_p))

    select case (os_p%variables_set)
    case ('GYRE')
       vr%set = GYRE_SET
    case ('DZIEM')
       vr%set = DZIEM_SET
    case ('JCD')
       vr%set = JCD_SET
    case ('MIX')
       vr%set = MIX_SET
    case ('LAGP')
       vr%set = LAGP_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    ! Finish

    return

  end function rad_vars_t_

  !****

  subroutine stencil (this, pt)

    class(rad_vars_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    call check_model(this%ml, [I_V_2,I_U,I_C_1])

    n_s = SIZE(pt)

    if (ALLOCATED(this%coeffs)) deallocate(this%coeffs)
    allocate(this%coeffs(n_s,J_LAST))

    do i = 1, n_s
       if (this%ml%is_vacuum(pt(i))) then
          $ASSERT(this%set /= LAGP_SET,Cannot use LAGP variables at vacuum points)
          this%coeffs(i,J_V_2) = HUGE(0._WP)
          this%coeffs(i,J_DV_2) = HUGE(0._WP)
       else
          this%coeffs(i,J_V_2) = this%ml%coeff(I_V_2, pt(i))
          this%coeffs(i,J_DV_2) = this%ml%dcoeff(I_V_2, pt(i))
       endif
       this%coeffs(i,J_C_1) = this%ml%coeff(I_C_1, pt(i))
       this%coeffs(i,J_DC_1) = this%ml%dcoeff(I_C_1, pt(i))
    end do

    ! Set up stencil for the rt component

    call this%rt%stencil(pt)

    ! Finish

    return

  end subroutine stencil

  !****

  function G (this, i, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%set)
    case (GYRE_SET)
       G = identity_matrix(2)
    case (DZIEM_SET)
       G = identity_matrix(2)
    case (JCD_SET)
       G = this%G_jcd_(i, omega)
    case (MIX_SET)
       G = identity_matrix(2)
    case (LAGP_SET)
       G = this%G_lagp_(i, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return
    
  end function G
  
  !****

  function G_jcd_ (this, i, omega) result (G)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    associate( &
         c_1 => this%coeffs(i,J_C_1))

      omega_c = this%rt%omega_c(i, omega)

      ! Set up the matrix

      G(1,1) = 1._WP
      G(1,2) = 0._WP

      G(2,1) = 0._WP
      G(2,2) = 1._WP/(c_1*omega_c**2)

    end associate
      
    ! Finish

    return

  end function G_jcd_

  !****

  function G_lagp_ (this, i, omega) result (G)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    associate( &
         V_2 => this%coeffs(i,J_V_2))

      ! Set up the matrix

      G(1,1) = 1._WP
      G(1,2) = 0._WP

      G(2,1) = -V_2
      G(2,2) = V_2

    end associate

    ! Finish

    return

  end function G_lagp_

  !****

  function H (this, i, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (GYRE_SET)
       H = identity_matrix(2)
    case (DZIEM_SET)
       H = identity_matrix(2)
    case (JCD_SET)
       H = this%H_jcd_(i, omega)
    case (MIX_SET)
       H = identity_matrix(2)
    case (LAGP_SET)
       H = this%H_lagp_(i, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function H

  !****

  function H_jcd_ (this, i, omega) result (H)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    associate( &
         c_1 => this%coeffs(i,J_C_1))

      omega_c = this%rt%omega_c(i, omega)

      ! Set up the matrix
      
      H(1,1) = 1._WP
      H(1,2) = 0._WP
      
      H(2,1) = 0._WP
      H(2,2) = c_1*omega_c**2

    end associate

    ! Finish

    return

  end function H_jcd_

  !****

  function H_lagp_ (this, i, omega) result (H)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    associate( &
         V_2 => this%coeffs(i,J_V_2))

      ! Set up the matrix

      H(1,1) = 1._WP
      H(1,2) = 0._WP
      
      H(2,1) = 1._WP
      H(2,2) = 1._WP/V_2

    end associate

    ! Finish

    return

  end function H_lagp_

  !****

  function dH (this, i, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case (GYRE_SET)
       dH = 0._WP
    case (DZIEM_SET)
       dH = 0._WP
    case (JCD_SET)
       dH = this%dH_jcd_(i, omega)
    case (MIX_SET)
       dH = 0._WP
    case (LAGP_SET)
       dH = this%dH_lagp_(i, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dH

  !****

  function dH_jcd_ (this, i, omega) result (dH)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    real(WP) :: omega_c

    ! Evaluate the derivative x dH/dx of the JCD-variables
    ! transformation matrix H

    associate( &
         c_1 => this%coeffs(i,J_C_1), &
         dc_1 => this%coeffs(i,J_DC_1))

      omega_c = this%rt%omega_c(i, omega)

      ! Set up the matrix (nb: the derivative of omega_c is neglected;
      ! this is incorrect when rotation is non-zero)
      
      dH(1,1) = 0._WP
      dH(1,2) = 0._WP

      dH(2,1) = 0._WP
      dH(2,2) = c_1*dc_1*omega_c**2

    end associate

    ! Finish

    return

  end function dH_jcd_

  !****

  function dH_lagp_ (this, i, omega) result (dH)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: i
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    ! Evaluate the derivative x dH/dx of the LAGP-variables
    ! transformation matrix H

    associate( &
         V_2 => this%coeffs(i,J_V_2), &
         dV_2 => this%coeffs(i,J_DV_2))

      ! Set up the matrix

      dH(1,1) = 0._WP
      dH(1,2) = 0._WP

      dH(2,1) = 0._WP
      dH(2,2) = -dV_2/V_2

    end associate

    ! Finish

    return

  end function dH_lagp_
  
end module gyre_rad_vars
