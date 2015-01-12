! Module   : gyre_cimplex
! Purpose  : complex root finding using downhill simplex minimization
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

module gyre_cimplex

  ! Uses

  use core_kinds

  use gyre_ext
  use gyre_extfunc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: cimplex_t
     private
     type(c_ext_t)                   :: cx(3)
     type(c_ext_t)                   :: f_cx(3)
     integer                         :: i_hi
     integer                         :: i_lo
     class(c_extfunc_t), allocatable :: cf
   contains
     private
     procedure, public :: refine => refine_
     procedure, public :: lowest => lowest_
  end type cimplex_t

  ! Interfaces

  interface cimplex_t
     module procedure cimplex_t_
     module procedure cimplex_t_verts_
  end interface cimplex_t

  ! Access specifiers

  private

  public :: cimplex_t

contains

  function cimplex_t_ (cf) result (cm)

    class(c_extfunc_t), intent(in) :: cf
    type(cimplex_t)                :: cm
       
    ! Construct the cimplex_t

    allocate(cm%cf, SOURCE=cf)

    ! Finish

    return

  end function cimplex_t_

!****

  function cimplex_t_verts_ (cx, f_cx, cf) result (cm)

    class(c_ext_t), intent(in)     :: cx(:)
    class(c_ext_t), intent(in)     :: f_cx(:)
    class(c_extfunc_t), intent(in) :: cf
    type(cimplex_t)                :: cm
       
    $CHECK_BOUNDS(SIZE(cx),3)
    $CHECK_BOUNDS(SIZE(f_cx),3)

    ! Construct the cimplex_t

    cm = cimplex_t(cf)

    ! Set up the starting vertices

    cm%cx = cx
    cm%f_cx = f_cx

    ! Find the highest and lowest vertices

    cm%i_hi = absmaxloc_(cm%f_cx)
    cm%i_lo = absminloc_(cm%f_cx)

    ! Finish

    return

  end function cimplex_t_verts_

!****

  subroutine refine_ (this, toler, n_iter)

    class(cimplex_t), intent(inout)  :: this
    real(WP), intent(in)             :: toler
    integer, optional, intent(inout) :: n_iter

    real(WP), parameter :: EPS = 4._WP*EPSILON(0._WP)
    logical, parameter  :: VERBOSE = .FALSE.

    integer         :: n_iter_
    real(WP)        :: tol
    type(cimplex_t) :: cm
    type(cimplex_t) :: cm_new
    integer         :: i
    integer         :: j

    if(PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 500
    end if

    ! Iterate the cimplex using the downhill simplex (Nelder-Mead)
    ! algorithm, until convergence or the number of iterations is
    ! exceeded. This implementation follows the same approach as
    ! described in Wikipedia

    tol = MAX(toler, EPS)

    ! ! This appears to be a gfortran bug; intrinsic assignment should work

    cm = cimplex_t(this%cf)

    cm%cx = this%cx
    cm%f_cx = this%f_cx

    cm%i_lo = this%i_lo
    cm%i_hi = this%i_hi

    if(VERBOSE) then

       write(OUTPUT_UNIT, *) 'Initial vertices:'
       do i = 1,3
          write(OUTPUT_UNIT, 100) i, fraction(this%cx(i)), exponent(this%cx(i)), fraction(this%f_cx(i)), exponent(this%f_cx(i))
100       format(I0,2(1X,E16.8),1X,I0,1X,E16.8,1X,I0)
       end do
    endif

    iterate_loop : do j = 1, n_iter_

       ! Check for convergence (equal minima)

       if (ABS(cm%f_cx(cm%i_hi))-ABS(cm%f_cx(cm%i_lo)) < &
            tol*ABS(cm%f_cx(cm%i_hi))) exit iterate_loop
          
       ! Update the cimplex

       ! Reflect

       cm_new = ooze_(cm, -1._WP, cm%i_hi)

       if (cm_new%i_lo == cm%i_hi) then

          ! Expand

          cm_new = ooze_(cm_new, 2._WP, cm%i_hi)

       elseif (cm_new%i_hi == cm%i_hi) then

          ! Contract

          cm_new = ooze_(cm, 0.5_WP, cm%i_hi)

          if (cm_new%i_hi == cm%i_hi) then

             ! Reduce

             cm_new = shrink_(cm, 0.5_WP, cm%i_lo)

          endif

       endif

       ! Check for convergence (change in coordinates)

       if (ABS(cm%cx(cm%i_hi) - cm%cx(cm%i_lo)) < &
            tol*MAX(ABS(cm%cx(cm%i_hi)), ABS(cm%cx(cm%i_lo)))) exit iterate_loop

       cm = cm_new

       if (VERBOSE) then
          write(OUTPUT_UNIT, *) 'Current vertices:'
          do i = 1,3
             write(OUTPUT_UNIT, *) i, cmplx(cm%cx(i)), real(cm%f_cx(i))
          end do
       endif

    end do iterate_loop

    if (PRESENT(n_iter)) then
       n_iter = j
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    if (VERBOSE) then
       write(OUTPUT_UNIT, *) 'Final vertices:'
       do i = 1,3
          write(OUTPUT_UNIT, *) i, cmplx(cm%cx(i)), real(cm%f_cx(i))
       end do
    endif

    ! this%cimplex_t = cm

    this%cx = cm%cx
    this%f_cx = cm%f_cx

    this%i_lo = cm%i_lo
    this%i_hi = cm%i_hi

    ! Finish

    return

  end subroutine refine_

!****

  function lowest_ (this) result (cx)

    class(cimplex_t), intent(in) :: this
    type(c_ext_t)                :: cx

    ! Return the lowest point of the cimplex

    cx = this%cx(this%i_lo)

    ! Finish

    return

  end function lowest_
  
!****

  function ooze_ (cm, scale, i) result (cm_new)

    type(cimplex_t), intent(inout) :: cm
    real(WP), intent(in)           :: scale
    integer, intent(in)            :: i
    type(cimplex_t)                :: cm_new

    integer       :: j
    type(c_ext_t) :: cx_cen
    type(c_ext_t) :: cx_new
    type(c_ext_t) :: f_cx_new

    ! (Provisionally) ooze the cimplex by scaling vertex i

    ! Find the centroid excluding vertex i

    cx_cen = c_ext_t(0._WP)

    do j = 1, 3
       if (j /= i) then
          cx_cen = cx_cen + cm%cx(j)
       endif
    enddo

    cx_cen = cx_cen/2._WP

    ! Calculate the new vertex position & value

    cx_new = cx_cen + scale*(cm%cx(i) - cx_cen)
    f_cx_new = cm%cf%eval(cx_new)

    ! Accept or reject the new position

    cm_new = cm

    if (ABS(f_cx_new) < ABS(cm%f_cx(i))) then

       cm_new%cx(i) = cx_new
       cm_new%f_cx(i) = f_cx_new

       cm_new%i_lo = absminloc_(cm_new%f_cx)
       cm_new%i_hi = absmaxloc_(cm_new%f_cx)

    endif

    ! Finish

    return

  end function ooze_

!****

  function shrink_ (cm, scale, i) result (cm_new)

    type(cimplex_t), intent(in) :: cm
    real(WP), intent(in)        :: scale
    integer, intent(in)         :: i
    type(cimplex_t)             :: cm_new

    integer :: j

    ! Shrink the cimplex by scale, keeping vertex i fixed

    cm_new = cm

    do j = 1,3
       if (j /= i) then
          cm_new%cx(j) = cm%cx(i) + scale*(cm%cx(j) - cm%cx(i))
          cm_new%f_cx(j) = cm_new%cf%eval(cm_new%cx(j))
       endif
    end do
          
    cm_new%i_lo = absminloc_(cm_new%f_cx)
    cm_new%i_hi = absmaxloc_(cm_new%f_cx)
    
    ! Finish

    return

  end function shrink_

!****

  function absmaxloc_ (cx) result (i_max)

    type(c_ext_t), intent(in) :: cx(:)
    integer                   :: i_max

    real(WP) :: f(SIZE(cx))
    integer  :: e(SIZE(cx))

    ! Return the index of the maximum absolute value in cx

    f = FRACTION(ABS(cx))
    e = MERGE(EXPONENT(ABS(cx)), -HUGE(0), f > 0._WP)

    i_max = MAXLOC(f, MASK=e==MAXVAL(e), DIM=1)

    ! Finish

    return

  end function absmaxloc_

!****

  function absminloc_ (cx) result (i_min)

    type(c_ext_t), intent(in) :: cx(:)
    integer                   :: i_min

    real(WP) :: f(SIZE(cx))
    integer  :: e(SIZE(cx))

    ! Return the index of the minimum absolute value in cx

    f = FRACTION(ABS(cx))
    e = MERGE(EXPONENT(ABS(cx)), -HUGE(0), f > 0._WP)

    i_min = MINLOC(f, MASK=e==MINVAL(e), DIM=1)

    ! Finish

    return

  end function absminloc_

end module gyre_cimplex
