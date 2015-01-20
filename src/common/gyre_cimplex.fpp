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
     integer                         :: i(3)
     class(c_extfunc_t), allocatable :: cf
   contains
     private
     procedure         :: set_indices => set_indices_
     procedure, public :: verts => verts_
     procedure, public :: values => values_
  end type cimplex_t

  ! Interfaces

  interface cimplex_t
     module procedure cimplex_t_
     module procedure cimplex_t_verts_
  end interface cimplex_t

  interface refine
     module procedure refine_
  end interface refine

  ! Access specifiers

  private

  public :: cimplex_t
  public :: refine

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

  function cimplex_t_verts_ (cf, cx, f_cx) result (cm)

    class(c_extfunc_t), intent(in)       :: cf
    class(c_ext_t), intent(in)           :: cx(:)
    class(c_ext_t), optional, intent(in) :: f_cx(:)
    type(cimplex_t)                      :: cm
       
    $CHECK_BOUNDS(SIZE(cx),3)
    $CHECK_BOUNDS(SIZE(f_cx),3)

    ! Construct the cimplex_t

    cm = cimplex_t(cf)

    ! Set up the starting vertices

    cm%cx = cx

    if (PRESENT(f_cx)) then
       cm%f_cx = f_cx
    else
       cm%f_cx(1) = cm%cf%eval(cx(1))
       cm%f_cx(2) = cm%cf%eval(cx(2))
       cm%f_cx(3) = cm%cf%eval(cx(3))
    end if

    call cm%set_indices()

    ! Finish

    return

  end function cimplex_t_verts_

!****

  subroutine set_indices_ (this)

    class(cimplex_t), intent(inout) :: this

    integer       :: i(3)
    type(r_ext_t) :: a(3)

    ! Set up the ranking indices (yes, this uses a bubble sort!)

    i = [1,2,3]

    a = ABS(this%f_cx)

    if (a(i(1)) > a(i(2))) call swap_(1, 2)
    if (a(i(2)) > a(i(3))) call swap_(2, 3)
    if (a(i(1)) > a(i(2))) call swap_(1, 2)

    this%i = i

    ! Finish

    return

  contains

    subroutine swap_ (j_a, j_b)

      integer, intent(in) :: j_a
      integer, intent(in) :: j_b

      integer :: i_tmp

      i_tmp = i(j_a)
      i(j_a) = i(j_b)
      i(j_b) = i_tmp

      return

    end subroutine swap_

  end subroutine set_indices_

!****

  function verts_ (this) result (cx)

    class(cimplex_t), intent(in) :: this
    type(c_ext_t)                :: cx(3)

    ! Get the vertices, sorted into increasing (value) order

    cx = this%cx(this%i)

    ! Finish

    return

  end function verts_

!****

  function values_ (this) result (f_cx)

    class(cimplex_t), intent(in) :: this
    type(c_ext_t)                :: f_cx(3)

    ! Get the values, sorted into increasing order

    f_cx = this%f_cx(this%i)

    ! Finish

    return

  end function values_

!****

  subroutine refine_ (cm, cx_tol, n_iter, relative_tol)

    type(cimplex_t), intent(inout)   :: cm
    type(r_ext_t), intent(in)        :: cx_tol
    integer, optional, intent(inout) :: n_iter
    logical, optional, intent(in)    :: relative_tol

    logical, parameter  :: VERBOSE = .FALSE.

    logical         :: relative_tol_
    integer         :: i_iter
    type(r_ext_t)   :: tol
    type(cimplex_t) :: cm_new
    integer         :: i

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Refine the cimplex until convergence or the number of iterations
    ! is exceeded. This implementation follows the Nelder-Mead
    ! algorithm as described in Wikipedia

    if (VERBOSE) then

       write(OUTPUT_UNIT, *) 'Initial vertices:'
       do i = 1,3
          write(OUTPUT_UNIT, 100) i, fraction(cm%cx(i)), exponent(cm%cx(i)), fraction(cm%f_cx(i)), exponent(cm%f_cx(i))
100       format(I0,2(1X,E16.8),1X,I0,1X,E16.8,1X,I0)
       end do
    endif

    i_iter = 0

    iterate_loop : do

       if (PRESENT(n_iter)) then
          if (i_iter >= n_iter) exit iterate_loop
       end if

       i_iter = i_iter + 1

       ! Check for convergence (equal minima)

       !if (ABS(cm%f_cx(cm%i(3)))-ABS(cm%f_cx(cm%i(1))) < &
       !     tol*ABS(cm%f_cx(cm%i(3)))) exit iterate_loop
          
       ! Update the cimplex

       ! Reflect

       cm_new = ooze_(cm, -1._WP, cm%i(3))

       if (cm_new%i(1) == cm%i(3)) then

          ! Expand

          cm_new = ooze_(cm_new, 2._WP, cm%i(3))

       elseif (cm_new%i(3) == cm%i(3)) then

          ! Contract

          cm_new = ooze_(cm, 0.5_WP, cm%i(3))

          if (cm_new%i(3) == cm%i(3)) then

             ! Reduce

             cm_new = shrink_(cm, 0.5_WP, cm%i(1))

          endif

       endif

       ! Check for convergence (change in coordinates)

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + cx_tol)*MAX(ABS(cm%cx(cm%i(3))), ABS(cm%cx(cm%i(1))))
       else
          tol = 4._WP*EPSILON(0._WP)*MAX(ABS(cm%cx(cm%i(3))), ABS(cm%cx(cm%i(1)))) + cx_tol
       endif

       if (ABS(cm%cx(cm%i(3)) - cm%cx(cm%i(1))) <= tol) exit iterate_loop

       cm = cm_new

       if (VERBOSE) then
          write(OUTPUT_UNIT, *) 'Current vertices:'
          do i = 1,3
             write(OUTPUT_UNIT, *) i, cmplx(cm%cx(i)), real(cm%f_cx(i))
          end do
       endif

    end do iterate_loop

    if (PRESENT(n_iter)) then
       n_iter = i_iter
    endif

    if (VERBOSE) then
       write(OUTPUT_UNIT, *) 'Final vertices:'
       do i = 1,3
          write(OUTPUT_UNIT, *) i, cmplx(cm%cx(i)), real(cm%f_cx(i))
       end do
    endif

    ! Finish

    return

  end subroutine refine_

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

       call cm_new%set_indices()

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

    call cm_new%set_indices()
    
    ! Finish

    return

  end function shrink_

end module gyre_cimplex
