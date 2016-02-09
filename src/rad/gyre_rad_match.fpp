! Module   : gyre_rad_match
! Purpose  : adiabatic radial match conditions
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

module gyre_rad_match

  ! Uses

  use core_kinds

  use gyre_rad_vars
  use gyre_diff
  use gyre_ext
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_diff_t) :: rad_match_t
     private
     class(model_t), pointer :: ml => null()
     type(rad_vars_t)        :: vr
     integer                 :: s
     real(WP)                :: x
   contains
     private
     procedure, public :: build
  end type rad_match_t

  ! Interfaces

  interface rad_match_t
     module procedure rad_match_t_
  end interface rad_match_t
  
  ! Access specifiers

  private
  public :: rad_match_t

contains

  function rad_match_t_ (ml, s_l, x_l, s_r, x_r, md_p, os_p) result (mt)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s_l
    real(WP)                            :: x_l
    integer, intent(in)                 :: s_r
    real(WP)                            :: x_r
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_match_t)                   :: mt

    $ASSERT(s_r == s_l+1,Invalid segment jump at match point)
    $ASSERT(x_r == x_l,Segments do not join at match point)

    ! Construct the rad_match_t

    mt%ml => ml

    mt%vr = rad_vars_t(ml, md_p, os_p)

    mt%s = s_l
    mt%x = x_l

    mt%n_e = 2

    ! Finish

    return

  end function rad_match_t_
    
  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(rad_match_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: E_l(:,:)
    real(WP), intent(out)          :: E_r(:,:)
    type(r_ext_t), intent(out)     :: scl

    real(WP) :: U_l
    real(WP) :: U_r

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Build the difference equations

    ! Calculate coefficients

    associate (s => this%s, &
               x => this%x)

      U_l = this%ml%U(s  , x)
      U_r = this%ml%U(s+1, x)

    end associate

    ! Evaluate the match conditions (y_1 continuous, y_2 not)

    E_l(1,1) = -1._WP
    E_l(1,2) = 0._WP
    
    E_l(2,1) = U_l
    E_l(2,2) = -U_l

    !

    E_r(1,1) = 1._WP
    E_r(1,2) = 0._WP

    E_r(2,1) = -U_r
    E_r(2,2) = U_r

    scl = r_ext_t(1._WP)

    ! Apply the variables transformation

    E_l = MATMUL(E_l, this%vr%H(this%s  , this%x, omega))
    E_r = MATMUL(E_r, this%vr%H(this%s+1, this%x, omega))

    ! Finish

    return

  end subroutine build

end module gyre_rad_match
