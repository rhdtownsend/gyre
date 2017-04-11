! Module   : gyre_diff_factory
! Purpose  : factory procedures for r_diff_t and c_diff_t types
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

module gyre_diff_factory

  ! Uses

  use core_kinds

  use gyre_colloc_diff
  use gyre_diff
  use gyre_eqns
  use gyre_magnus_diff
  use gyre_mirk_diff
  use gyre_num_par
  use gyre_point
  use gyre_trapz_diff

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_diff_t
     module procedure r_diff_t_
  end interface r_diff_t

  interface c_diff_t
     module procedure c_diff_t_
  end interface c_diff_t

  ! Access specifiers

  private

  public :: r_diff_t
  public :: c_diff_t

  ! Procedures

contains

  $define $DIFF_T $sub

  $local $T $1

  function ${T}_diff_t_ (eq, pt_a, pt_b, nm_p) result (df)

    class(${T}_eqns_t), intent(in)  :: eq
    type(point_t), intent(in)       :: pt_a
    type(point_t), intent(in)       :: pt_b
    type(num_par_t), intent(in)     :: nm_p
    class(${T}_diff_t), allocatable :: df
    
    ! Construct the ${T}_diff_t

    select case (nm_p%diff_scheme)
    case ('MAGNUS_GL2')
       allocate(df, SOURCE=${T}_magnus_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(df, SOURCE=${T}_magnus_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(df, SOURCE=${T}_magnus_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('COLLOC_GL2')
       allocate(df, SOURCE=${T}_colloc_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('COLLOC_GL4')
       allocate(df, SOURCE=${T}_colloc_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('COLLOC_GL6')
       allocate(df, SOURCE=${T}_colloc_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('TRAPZ')
       allocate(df, SOURCE=${T}_trapz_diff_t(eq, pt_a, pt_b, SPREAD(0.5_WP, DIM=1, NCOPIES=eq%n_e)))
    case ('MIRK')
       allocate(df, SOURCE=${T}_mirk_diff_t(eq, pt_a, pt_b))
    case default
       $ABORT(Invalid diff_scheme)
    end select

    ! Finish

    return

  end function ${T}_diff_t_

  $endsub

  $DIFF_T(r)
  $DIFF_T(c)

end module gyre_diff_factory
