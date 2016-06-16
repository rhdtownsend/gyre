! Incfile  : gyre_cheby
! Purpose  : Chebyshev interpolation & expansions
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_cheby

  ! Uses

  use core_kinds
  use core_constants
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_memory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: cheby_t
!     private
     real(WP), allocatable :: f(:)
     real(WP), allocatable :: c(:)
     real(WP)              :: x_a
     real(WP)              :: x_b
     integer               :: n
   contains 
     private
     procedure       :: eval_r_
     procedure       :: eval_c_
     generic, public :: eval => eval_r_, eval_c_
  end type cheby_t

  ! Interfaces

  interface cheby_t
     module procedure cheby_t_func_
     module procedure cheby_t_coeffs_
     module procedure cheby_t_tol_
  end interface cheby_t
  
  $if ($HDF5)
  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write
  $endif

  ! Access specifiers

  private

  public :: cheby_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif

  ! Procedures

contains

  function cheby_t_func_ (x_a, x_b, n, func) result (cb)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    integer, intent(in)  :: n
    interface
       function func (x)
         use core_kinds
         real(WP), intent(in) :: x
         real(WP)             :: func
       end function func
    end interface
    type(cheby_t)        :: cb

    real(WP) :: f(n+1)
    integer  :: j
    real(WP) :: x
    real(WP) :: u

    ! Construct the cheby_t of degree n, by sampling the function at
    ! the n+1 extremal points of T_n

    do j = 1, n+1

       if (j == 1) then
          x = x_a
       elseif (j == n+1) then
          x = x_b
       else
          u = COS((j-1)*PI/n)
          x = 0.5_WP*((1._WP+u)*x_a + (1._WP-u)*x_b)
       endif

       f(j) = func(x)

    end do

    cb = cheby_t_vals_(x_a, x_b, f)

    ! Finish

    return

  end function cheby_t_func_

  !****

  function cheby_t_tol_ (x_a, x_b, tol, func) result (cb)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP), intent(in) :: tol
    interface
       function func (x)
         use core_kinds
         real(WP), intent(in) :: x
         real(WP)             :: func
       end function func
    end interface
    type(cheby_t)        :: cb

    integer, parameter  :: N_0 = 16
    integer, parameter  :: M = 8
    real(WP), parameter :: EPS = EPSILON(1._WP)

    integer  :: n
    integer  :: j
    real(WP) :: toler

    ! Construct a cheby_t by choosing an n such that all high-order
    ! (neglected) coefficients are below a (relative) tolerance toler

    ! Starting n

    n = N_0

    ! Increase n until at least M trailing coefficients of the cheby_t
    ! are below toler

    do

       cb = cheby_t_func_(x_a, x_b, n, func)

       toler = (tol + 10._WP*EPS)*MAXVAL(ABS(cb%f))

       do j = n+1, 1, -1
          if (cb%c(j) > toler) exit
       end do

       if (n+1-j >= M) exit

       n = n*2

    end do

    ! Re-create the cheby_t with the optimal n

    cb = cheby_t_func_(x_a, x_b, j, func)

    ! Finish

    return

  end function cheby_t_tol_

  !****

  function cheby_t_vals_ (x_a, x_b, f) result (cb)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP), intent(in) :: f(:)
    type(cheby_t)        :: cb

    ! Construct the cheby_t of degree n, using the supplied function
    ! values at the extremal points of T_n

    cb%f = f
    cb%c = c_from_f(f)

    cb%x_a = x_a
    cb%x_b = x_b
    
    cb%n = SIZE(f) - 1
    
    ! Finish

    return

  end function cheby_t_vals_

  !****

  function cheby_t_coeffs_ (x_a, x_b, c) result (cb)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP), intent(in) :: c(:)
    type(cheby_t)        :: cb

    ! Construct the cheby_t of degree n, using the supplied expansion
    ! coefficients

    cb%c = c
    cb%f = f_from_c(c)

    cb%x_a = x_a
    cb%x_b = x_b
    
    cb%n = SIZE(c)

    ! Finish

    return

  end function cheby_t_coeffs_

  !****

  $if ($HDF5)

  subroutine read_ (hg, cb)

    type(hgroup_t), intent(inout) :: hg
    type(cheby_t), intent(out)    :: cb

    real(WP)              :: x_a
    real(WP)              :: x_b
    real(WP), allocatable :: f(:)

    ! Read the cheby_t

    call read_attr(hg, 'x_a', x_a)
    call read_attr(hg, 'x_b', x_b)

    call read_dset_alloc(hg, 'f', f)

    cb = cheby_t_vals_(x_a, x_b, f)

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, cb)

    type(hgroup_t), intent(inout) :: hg
    type(cheby_t), intent(in)     :: cb

    ! Write the cheby_t

    call write_attr(hg, 'x_a', cb%x_a)
    call write_attr(hg, 'x_b', cb%x_b)

    call write_dset(hg, 'f', cb%f)

    ! Finish

    return

  end subroutine write_

  $endif

  !****

  $define $EVAL $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function eval_${SUFFIX}_ (this, x) result (f)

    class(cheby_t), intent(in) :: this
    $TYPE(WP), intent(in)      :: x
    $TYPE(WP)                  :: f

    $TYPE(WP) :: u
    $TYPE(WP) :: l
    $TYPE(WP) :: s
    integer   :: j
    real(WP)  :: u_j
    real(WP)  :: w

    ! Evaluate the cheby_t at x, using Barycentric interpolation
    ! following eqn. 5.9 of Trefethen (Approximation Theory &
    ! Approximation Practice)

    if (x == this%x_a) then
       u = 1._WP
    elseif (x == this%x_b) then
       u = -1._WP
    else
       u = (2._WP*x - (this%x_a + this%x_b))/(this%x_a - this%x_b)
    endif

    l = 1._WP
    s = 0._WP

    do j = 1, this%n+1

       if (j == 1) then
          u_j = 1._WP
       elseif (j == this%n+1) then
          u_j = -1._WP
       else
          u_j = COS((j-1)*PI/this%n)
       endif

       if (u == u_j) then
          f = this%f(j)
          return
       endif

       ! Calculate weights w =lambda / [2**(n-1)/n]

       if (j == 1 .OR. j == this%n+1) then
          w = 0.5_WP*(-1._WP)**(j-1)
       else
          w = (-1._WP)**(j-1)
       endif

       ! Update sum s and product l
       
       s = s + w*this%f(j)/(u - u_j)
       l = l*(u - u_j)

    end do

    f = l*s*(2._WP)**(this%n-1)/this%n

    ! Finish

    return

  end function eval_${SUFFIX}_

  $endsub

  $EVAL(r,real)
  $EVAL(c,complex)
  
  function c_from_f (f) result (c)

    real(WP), intent(in) :: f(:)
    real(WP)             :: c(SIZE(f))

    integer  :: n
    integer  :: k
    integer  :: j
    real(WP) :: v
    real(WP) :: w

    ! Calculate the Chebyshev expansion coefficients c from the
    ! function f sampled at the extremal points of T_n

    n = SIZE(c) - 1
   
    do k = 1, n+1

       c(k) = 0._WP

       do j = 1, n+1

          if (j == 1) then
             v = 1._WP
             w = 0.5_WP
          elseif (j == n+1) then
             v = (-1._WP)**(k-1)
             w = 0.5_WP
          else
             v = COS((j-1)*(k-1)*PI/n)
             w = 1._WP
          endif
          
          c(k) = c(k) + w*v*f(j)

       end do

       if (k == 1 .OR. k == n+1) then
          c(k) = c(k)/n
       else
          c(k) = 2._WP*c(k)/n
       end if

    end do

    ! Finish

    return

  end function c_from_f

  !****

  function f_from_c (c) result (f)

    real(WP), intent(in) :: c(:)
    real(WP)             :: f(SIZE(c))

    integer  :: n
    integer  :: j
    integer  :: k
    real(WP) :: v

    ! Calculate the function f sampled at the extremal points of T_n,
    ! from the Chebyshev expansion coefficients c

    n = SIZE(c)-1

    do j = 1, n+1

       f(j) = 0._WP

       do k = 1, n+1

          if (k == 1) then
             v = 1._WP
          elseif (k == n+1) then
             v = (-1._WP)**(j-1)
          else
             v = COS((j-1)*(k-1)*PI/n)
          endif

          f(j) = f(j) + v*c(k)

       end do

    end do

    ! Finish

    return

  end function f_from_c

end module gyre_cheby
