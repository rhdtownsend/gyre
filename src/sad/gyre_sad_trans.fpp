! Module   : gyre_sad_vars
! Purpose  : static adiabatic variables transformations
!
! Copyright 2019 Rich Townsend
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

module gyre_sad_trans

  ! Uses

  use core_kinds

  use gyre_context
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_state
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GYRE_SET = 0

  ! Derived-type definitions

  type :: sad_trans_t
     private
     type(context_t), pointer :: cx => null()
     integer                  :: set
     integer                  :: n_e
   contains
     private
     procedure, public :: stencil
     procedure, public :: trans_eqns
     procedure, public :: trans_cond
     procedure, public :: trans_vars
     procedure         :: G_
     procedure         :: H_
     procedure         :: dH_
  end type sad_trans_t

  ! Interfaces

  interface sad_trans_t
     module procedure sad_trans_t_
  end interface sad_trans_t

  ! Access specifiers

  private

  public :: sad_trans_t

  ! Procedures

contains

  function sad_trans_t_ (cx, md_p, os_p) result (tr)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(sad_trans_t)                    :: tr

    ! Construct the sad_trans_t

    tr%cx => cx

    select case (os_p%variables_set)
    case ('GYRE')
       tr%set = GYRE_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    tr%n_e = 2

    ! Finish

    return

  end function sad_trans_t_

  !****

  subroutine stencil (this, pt)

    class(sad_trans_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt(:)

    ! Finish

    return

  end subroutine stencil

  !****

  subroutine trans_eqns (this, xA, i, st, from)

    class(sad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: xA(:,:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)
    real(WP) :: dH(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(xA, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(xA, 2),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Transform equations to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET) then
          G = this%G_(i, st)
          H = this%H_(i, st)
          dH = this%dH_(i, st)
          xA = MATMUL(G, MATMUL(xA, H) - dH)
       endif

    else

       ! Convert to

       $ABORT(Not currently supported)

    end if

    ! Finish

    return

  end subroutine trans_eqns
  
  !****

  subroutine trans_cond (this, C, i, st, from)

    class(sad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: C(:,:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(C, 2),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Transform boundary/match conditions to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET) then
          H = this%H_(i, st)
          C = MATMUL(C, H)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET) then
          G = this%G_(i, st)
          C = MATMUL(C, G)
       endif

    end if

    ! Finish

    return

  end subroutine trans_cond

  !****

  subroutine trans_vars (this, y, i, st, from)

    class(sad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: y(:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(y),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Convert variables to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET) then
          G = this%G_(i, st)
          y = MATMUL(G, y)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET) then
          H = this%H_(i, st)
          y = MATMUL(H, y)
       endif

    end if

    ! Finish

    return

  end subroutine trans_vars

  !****

  function G_ (this, i, st) result (G)

    class(sad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables from
    ! GYRE's canonical form

    select case (this%set)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return
    
  end function G_
  
  !****

  function H_ (this, i, st) result (H)

    class(sad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function H_

  !****

  function dH_ (this, i, st) result (dH)

    class(sad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: dH(this%n_e,this%n_e)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dH_

end module gyre_sad_trans
