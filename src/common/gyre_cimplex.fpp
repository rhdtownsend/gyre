! Module   : gyre_cimplex
! Purpose  : complex root finding using downhill simplex minimization
!
! Copyright 2013 Rich Townsend
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

  use gyre_ext_arith
  use gyre_ext_func

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: cimplex_t
     type(ext_complex_t)            :: ez(3)
     type(ext_real_t)               :: f_ez(3)
     integer                        :: i_hi
     integer                        :: i_lo
     class(ext_func_t), allocatable :: ef
   contains
     private
     procedure, public :: refine => refine_
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

  function cimplex_t_ (ef) result (cm)

    class(ext_func_t), intent(in) :: ef
    type(cimplex_t)               :: cm
       
    ! Construct the cimplex_t

    allocate(cm%ef, SOURCE=ef)

    ! Finish

    return

  end function cimplex_t_

!****

  function cimplex_t_verts_ (ez, f_ez, ef) result (cm)

    class(ext_complex_t), intent(in) :: ez(:)
    class(ext_complex_t), intent(in) :: f_ez(:)
    class(ext_func_t), intent(in)    :: ef
    type(cimplex_t)                  :: cm
       
    integer  :: i

    $CHECK_BOUNDS(SIZE(ez),3)
    $CHECK_BOUNDS(SIZE(f_ez),3)

    ! Construct the cimplex_t

    cm = cimplex_t(ef)

    ! Set up the starting vertices

    cm%ez = ez
    cm%f_ez = ABS(f_ez)

    ! Find the highest and lowest vertices

    cm%i_hi = maxloc_pos_(cm%f_ez)
    cm%i_lo = minloc_pos_(cm%f_ez)

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
    cm = cimplex_t(this%ef)
    cm%ez = this%ez
    cm%f_ez = this%f_ez
    cm%i_lo = this%i_lo
    cm%i_hi = this%i_hi

    if(VERBOSE) then
       write(OUTPUT_UNIT, *) 'Initial vertices:'
       do i = 1,3
          write(OUTPUT_UNIT, 100) i, fraction(this%ez(i)), exponent(this%ez(i)), fraction(this%f_ez(i)), exponent(this%f_ez(i))
100       format(I0,2(1X,E16.8),1X,I0,1X,E16.8,1X,I0)
       end do
    endif

    iterate_loop : do j = 1, n_iter_

       ! Check for convergence (equal minima)

       if (cm%f_ez(cm%i_hi)-cm%f_ez(cm%i_lo) < &
            tol*cm%f_ez(cm%i_hi)) exit iterate_loop
          
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

       if (ABS(cm%ez(cm%i_hi) - cm%ez(cm%i_lo)) < &
            tol*MAX(ABS(cm%ez(cm%i_hi)), ABS(cm%ez(cm%i_lo)))) exit iterate_loop

       cm = cm_new

       if (VERBOSE) then
          write(OUTPUT_UNIT, *) 'Current vertices:'
          do i = 1,3
             write(OUTPUT_UNIT, *) i, cmplx(cm%ez(i)), real(cm%f_ez(i))
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
          write(OUTPUT_UNIT, *) i, cmplx(cm%ez(i)), real(cm%f_ez(i))
       end do
    endif

    this%ez = cm%ez
    this%f_ez = cm%f_ez

    this%i_lo = cm%i_lo
    this%i_hi = cm%i_hi

    ! Finish

    return

  end subroutine refine_

!****

  function ooze_ (cm, scale, i) result (cm_new)

    type(cimplex_t), intent(inout) :: cm
    real(WP), intent(in)           :: scale
    integer, intent(in)            :: i
    type(cimplex_t)                :: cm_new

    integer             :: j
    type(ext_complex_t) :: ez_cen
    type(ext_complex_t) :: ez_new
    type(ext_real_t)    :: f_ez_new

    ! (Provisionally) ooze the cimplex by scaling vertex i

    ! Find the centroid excluding vertex i

    ez_cen = ext_complex_t(0._WP)

    do j = 1, 3
       if (j /= i) then
          ez_cen = ez_cen + cm%ez(j)
       endif
    enddo

    ez_cen = ez_cen/2._WP

    ! Calculate the new vertex position & value

    ez_new = ez_cen + scale*(cm%ez(i) - ez_cen)
    f_ez_new = ABS(cm%ef%eval(ez_new))

    ! Accept or reject the new position

    cm_new = cm

    if (f_ez_new < cm%f_ez(i)) then

       cm_new%ez(i) = ez_new
       cm_new%f_ez(i) = f_ez_new

       cm_new%i_lo = minloc_pos_(cm_new%f_ez)
       cm_new%i_hi = maxloc_pos_(cm_new%f_ez)

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
          cm_new%ez(j) = cm%ez(i) + scale*(cm%ez(j) - cm%ez(i))
          cm_new%f_ez(j) = ABS(cm_new%ef%eval(cm_new%ez(j)))
       endif
    end do
          
    cm_new%i_lo = minloc_pos_(cm_new%f_ez)
    cm_new%i_hi = maxloc_pos_(cm_new%f_ez)
    
    ! Finish

    return

  end function shrink_

!****

  function maxloc_pos_ (ez) result (i_max)

    type(ext_real_t), intent(in) :: ez(:)
    integer                      :: i_max

    real(WP) :: f(SIZE(ez))
    integer  :: e(SIZE(ez))

    ! Return the index of the maximum value in ez, *assuming
    ! all-positive arguments*

    f = FRACTION(ez)
    e = MERGE(EXPONENT(ez), -HUGE(0), f > 0._WP)

    i_max = MAXLOC(f, MASK=e==MAXVAL(e), DIM=1)

    ! Finish

    return

  end function maxloc_pos_

!****

  function minloc_pos_ (ez) result (i_min)

    type(ext_real_t), intent(in) :: ez(:)
    integer                      :: i_min

    real(WP) :: f(SIZE(ez))
    integer  :: e(SIZE(ez))

    ! Return the index of the minimum value in ez, *assuming
    ! all-positive arguments*

    f = FRACTION(ez)
    e = MERGE(EXPONENT(ez), -HUGE(0), f > 0._WP)

    i_min = MINLOC(f, MASK=e==MINVAL(e), DIM=1)

    ! Finish

    return

  end function minloc_pos_

end module gyre_cimplex
