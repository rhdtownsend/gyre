! Module   : gyre_discrim
! Purpose  : discriminant evaluation
!
! Copyright 2021 Rich Townsend & The GYRE Team
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

module gyre_discrim

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface eval_discrim
     module procedure eval_discrim_r_
     module procedure eval_discrim_c_
  end interface eval_discrim

  ! Access specifiers

  private

  public :: eval_discrim

  ! Procedures

contains

  $define $DISCRIM $sub

  $local $T $1

  subroutine eval_discrim_${T}_ (bp, st, discrim)

    class(${T}_bvp_t), intent(inout) :: bp
    class(${T}_state_t), intent(in)  :: st
    type(${T}_ext_t)                 :: discrim

    ! Evaluate the discriminant

    call bp%build(st)
    call bp%factor()

    discrim = bp%det()

    ! Finish

    return

  end subroutine eval_discrim_${T}_

  $endsub

  $DISCRIM(r)
  $DISCRIM(c)

end module gyre_discrim
