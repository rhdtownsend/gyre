! Module   : gyre_nad_match
! Purpose  : nonadiabatic match conditions
!
! Copyright 2016 Rich Townsend
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

module gyre_nad_match

  ! Uses

  use core_kinds

  use gyre_diff
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_nad_vars
  use gyre_osc_par
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_diff_t) :: nad_match_t
     private
     class(model_t), pointer :: ml => null()
     type(nad_vars_t)        :: vr
     type(point_t)           :: pt_a
     type(point_t)           :: pt_b
   contains
     private
     procedure, public :: build
  end type nad_match_t

  ! Interfaces

  interface nad_match_t
     module procedure nad_match_t_
  end interface nad_match_t
  
  ! Access specifiers

  private
  public :: nad_match_t

contains

  function nad_match_t_ (ml, gr, k, md_p, os_p) result (mt)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    integer, intent(in)                 :: k
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_match_t)                   :: mt

    $ASSERT_DEBUG(k >= 1,Invalid index)
    $ASSERT_DEBUG(k < gr%n_k,Invalid index)

    $ASSERT_DEBUG(gr%pt(k+1)%s == gr%pt(k)%s+1,Mismatched segments)
    $ASSERT_DEBUG(gr%pt(k+1)%x == gr%pt(k)%x,Mismatched abscissae)

    ! Construct the nad_match_t

    mt%ml => ml

    mt%vr = nad_vars_t(ml, gr, md_p, os_p)

    mt%pt_a = gr%pt(k)
    mt%pt_b = gr%pt(k+1)

    mt%n_e = 6

    ! Finish

    return

  end function nad_match_t_
    
  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(nad_match_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: E_l(:,:)
    complex(WP), intent(out)       :: E_r(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP) :: V_l
    real(WP) :: V_r
    real(WP) :: U_l
    real(WP) :: U_r
    real(WP) :: nabla_ad_l
    real(WP) :: nabla_ad_r

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Build the difference equations

    ! Calculate coefficients

    associate (pt_a => this%pt_a, &
               pt_b => this%pt_b)

      V_l = this%ml%V_2(pt_a)*pt_a%x**2
      V_r = this%ml%V_2(pt_b)*pt_b%x**2

      U_l = this%ml%U(pt_a)
      U_r = this%ml%U(pt_b)

      nabla_ad_l = this%ml%nabla_ad(pt_a)
      nabla_ad_r = this%ml%nabla_ad(pt_b)

    end associate

    ! Evaluate the match conditions (y_1, y_3, y_6 continuous, y_2,
    ! y_4, y_5 not)

    E_l(1,1) = -1._WP
    E_l(1,2) = 0._WP
    E_l(1,3) = 0._WP
    E_l(1,4) = 0._WP
    E_l(1,5) = 0._WP
    E_l(1,6) = 0._WP
    
    E_l(2,1) = U_l
    E_l(2,2) = -U_l
    E_l(2,3) = 0._WP
    E_l(2,4) = 0._WP
    E_l(2,5) = 0._WP
    E_l(2,6) = 0._WP

    E_l(3,1) = 0._WP
    E_l(3,2) = 0._WP
    E_l(3,3) = -1._WP
    E_l(3,4) = 0._WP
    E_l(3,5) = 0._WP
    E_l(3,6) = 0._WP

    E_l(4,1) = -U_l
    E_l(4,2) = 0._WP
    E_l(4,3) = 0._WP
    E_l(4,4) = -1._WP
    E_l(4,5) = 0._WP
    E_l(4,6) = 0._WP

    E_l(5,1) = V_l*nabla_ad_l
    E_l(5,2) = -V_l*nabla_ad_l
    E_l(5,3) = 0._WP
    E_l(5,4) = 0._WP
    E_l(5,5) = -1._WP
    E_l(5,6) = 0._WP

    E_l(6,1) = 0._WP
    E_l(6,2) = 0._WP
    E_l(6,3) = 0._WP
    E_l(6,4) = 0._WP
    E_l(6,5) = 0._WP
    E_l(6,6) = -1._WP

    !

    E_r(1,1) = 1._WP
    E_r(1,2) = 0._WP
    E_r(1,3) = 0._WP
    E_r(1,4) = 0._WP
    E_r(1,5) = 0._WP
    E_r(1,6) = 0._WP

    E_r(2,1) = -U_r
    E_r(2,2) = U_r
    E_r(2,3) = 0._WP
    E_r(2,4) = 0._WP
    E_r(2,5) = 0._WP
    E_r(2,6) = 0._WP

    E_r(3,1) = 0._WP
    E_r(3,2) = 0._WP
    E_r(3,3) = 1._WP
    E_r(3,4) = 0._WP
    E_r(3,5) = 0._WP
    E_r(3,6) = 0._WP

    E_r(4,1) = U_r
    E_r(4,2) = 0._WP
    E_r(4,3) = 0._WP
    E_r(4,4) = 1._WP
    E_r(4,5) = 0._WP
    E_r(4,6) = 0._WP

    E_r(5,1) = -V_r*nabla_ad_r
    E_r(5,2) = V_r*nabla_ad_r
    E_r(5,3) = 0._WP
    E_r(5,4) = 0._WP
    E_r(5,5) = 1._WP
    E_r(5,6) = 0._WP

    E_r(6,1) = 0._WP
    E_r(6,2) = 0._WP
    E_r(6,3) = 0._WP
    E_r(6,4) = 0._WP
    E_r(6,5) = 0._WP
    E_r(6,6) = 1._WP

    scl = c_ext_t(1._WP)

    ! Apply the variables transformation

    E_l = MATMUL(E_l, this%vr%H(this%pt_a, omega))
    E_r = MATMUL(E_r, this%vr%H(this%pt_b, omega))

    ! Finish

    return

  end subroutine build

end module gyre_nad_match
