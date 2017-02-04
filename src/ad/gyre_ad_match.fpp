! Module   : gyre_ad_match
! Purpose  : adiabatic match conditions
!
! Copyright 2015-2017 Rich Townsend
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

module gyre_ad_match

  ! Uses

  use core_kinds

  use gyre_ad_vars
  use gyre_diff
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_U = 1
  
  integer, parameter :: J_LAST = J_U

  ! Derived-type definitions

  type, extends (r_diff_t) :: ad_match_t
     private
     class(model_t), pointer :: ml => null()
     type(ad_vars_t)         :: vr
     real(WP), allocatable   :: coeffs(:,:)
   contains
     private
     procedure         :: stencil_
     procedure, public :: build
  end type ad_match_t

  ! Interfaces

  interface ad_match_t
     module procedure ad_match_t_
  end interface ad_match_t
  
  ! Access specifiers

  private
  public :: ad_match_t

contains

  function ad_match_t_ (ml, pt_i, pt_a, pt_b, md_p, os_p) result (mt)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_a
    type(point_t), intent(in)           :: pt_b
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_match_t)                    :: mt

    $ASSERT_DEBUG(pt_a%s+1 == pt_b%s,Mismatched segments)
    $ASSERT_DEBUG(pt_a%x == pt_b%x,Mismatched abscissae)

    ! Construct the ad_match_t

    mt%ml => ml

    mt%vr = ad_vars_t(ml, pt_i, md_p, os_p)

    call mt%stencil_(pt_a, pt_b)

    mt%n_e = 4

    ! Finish

    return

  end function ad_match_t_
    
  !****

  subroutine stencil_ (this, pt_a, pt_b)

    class(ad_match_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt_a
    type(point_t), intent(in)        :: pt_b

    ! Calculate coefficients at the stencil points

    call check_model(this%ml, [I_U])

    allocate(this%coeffs(2,J_LAST))

    this%coeffs(1,J_U) = this%ml%coeff(I_U, pt_a)
    this%coeffs(2,J_U) = this%ml%coeff(I_U, pt_b)

    ! Set up stencil for the vr component

    call this%vr%stencil([pt_a,pt_b])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(ad_match_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP), intent(out)         :: E_l(:,:)
    real(WP), intent(out)         :: E_r(:,:)
    type(r_ext_t), intent(out)    :: scl

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Build the difference equations

    ! Calculate coefficients

    associate( &
      U_l => this%coeffs(1,J_U), &
      U_r => this%coeffs(2,J_U))

      ! Evaluate the match conditions (y_1, y_3 continuous, y_2, y_4
      ! not)

      E_l(1,1) = -1._WP
      E_l(1,2) = 0._WP
      E_l(1,3) = 0._WP
      E_l(1,4) = 0._WP
    
      E_l(2,1) = U_l
      E_l(2,2) = -U_l
      E_l(2,3) = 0._WP
      E_l(2,4) = 0._WP

      E_l(3,1) = 0._WP
      E_l(3,2) = 0._WP
      E_l(3,3) = -1._WP
      E_l(3,4) = 0._WP

      E_l(4,1) = -U_l
      E_l(4,2) = 0._WP
      E_l(4,3) = 0._WP
      E_l(4,4) = -1._WP

      !

      E_r(1,1) = 1._WP
      E_r(1,2) = 0._WP
      E_r(1,3) = 0._WP
      E_r(1,4) = 0._WP

      E_r(2,1) = -U_r
      E_r(2,2) = U_r
      E_r(2,3) = 0._WP
      E_r(2,4) = 0._WP
      
      E_r(3,1) = 0._WP
      E_r(3,2) = 0._WP
      E_r(3,3) = 1._WP
      E_r(3,4) = 0._WP

      E_r(4,1) = U_r
      E_r(4,2) = 0._WP
      E_r(4,3) = 0._WP
      E_r(4,4) = 1._WP

      scl = r_ext_t(1._WP)

    end associate

    ! Apply the variables transformation

    E_l = MATMUL(E_l, this%vr%H(1, omega))
    E_r = MATMUL(E_r, this%vr%H(2, omega))

    ! Finish

    return

  end subroutine build

end module gyre_ad_match
