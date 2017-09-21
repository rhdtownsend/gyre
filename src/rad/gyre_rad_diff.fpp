! Incfile  : gyre_rad_diff
! Purpose  : adiabatic radial difference equations
!
! Copyright 2016-2017 Rich Townsend
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

module gyre_rad_diff

  ! Uses

  use core_kinds

  use gyre_context
  use gyre_diff_factory
  use gyre_ext
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_eqns
  use gyre_rad_match
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_diff_t) :: rad_diff_t
     private
     class(r_diff_t), allocatable :: df
   contains
     private
     procedure, public :: build
  end type rad_diff_t

  ! Interfaces

  interface rad_diff_t
     module procedure rad_diff_t_
  end interface rad_diff_t

  ! Access specifiers

  private

  public :: rad_diff_t

  ! Procedures

contains

  function rad_diff_t_ (cx, pt_i, pt_a, pt_b, md_p, nm_p, os_p) result (df)

    type(context_t), pointer, intent(in) :: cx
    type(point_t), intent(in)            :: pt_i
    type(point_t), intent(in)            :: pt_a
    type(point_t), intent(in)            :: pt_b
    type(mode_par_t), intent(in)         :: md_p
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_diff_t)                     :: df

    type(rad_eqns_t) :: eq

    ! Construct the rad_diff_t

    if (pt_a%s == pt_b%s) then

       ! Regular subinterval; use difference equations

       eq = rad_eqns_t(cx, pt_i, md_p, os_p)
       
       allocate(df%df, SOURCE=r_diff_t(eq, pt_a, pt_b, nm_p))

    else

       ! Segment boundary; use match conditions

       allocate(df%df, SOURCE=rad_match_t(cx, pt_i, pt_a, pt_b, md_p, os_p))

    endif

    df%n_e = df%df%n_e

    ! Finish

    return

  end function rad_diff_t_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(rad_diff_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: E_l(:,:)
    real(WP), intent(out)         :: E_r(:,:)
    type(r_ext_t), intent(out)    :: scl

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Build the difference equations

    call this%df%build(st, E_l, E_r, scl)

    ! Finish

    return

  end subroutine build

end module gyre_rad_diff
