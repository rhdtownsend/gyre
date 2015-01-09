! Incfile  : gyre_cheby
! Purpose  : Chebyshev fitting functions
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

module gyre_cheby

  ! Uses

  use core_kinds
  use core_constants
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_memory

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: cheby_t
     private
     real(WP), allocatable :: c(:)
     real(WP)              :: x_a
     real(WP)              :: x_b
     integer               :: n
   contains 
     private
     procedure, public :: truncate => truncate_
     procedure         :: eval_r_
     procedure         :: eval_c_
     generic, public   :: eval => eval_r_, eval_c_
  end type cheby_t

  ! Interfaces

  interface cheby_t
     module procedure cheby_t_coeff_
     module procedure cheby_t_func_
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

  function cheby_t_coeff_ (x_a, x_b, c) result (cb)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP), intent(in) :: c(:)
    type(cheby_t)        :: cb

    ! Construct the cheby_t using the coefficients

    cb%c = c

    cb%x_a = x_a
    cb%x_b = x_b

    cb%n = SIZE(c)

    ! Finish

    return

  end function cheby_t_coeff_

!****

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

    real(WP) :: b_m
    real(WP) :: b_p
    integer  :: j
    integer  :: k
    real(WP) :: y
    real(WP) :: f(n)
    real(WP) :: fac
    real(DP) :: sum

    ! Construct the cheby_t, with n terms, from the function func
    ! defined on the domain [x_a,x_b]

    ! Set up the coefficients

    allocate(cb%c(n))

    b_m = 0.5_WP*(x_b - x_a)
    b_p = 0.5_WP*(x_b + x_a)

    do k = 1, n
       y = COS(PI*(k-0.5_DP)/n)
       f(k) = func(y*b_m + b_p)
    end do

    fac = 2._WP/n

    do j = 1, n
       sum = 0._DP
       do k = 1, n
          sum = sum + f(k)*COS((PI*(j-1))*((k-0.5_DP)/n))
       end do
       cb%c(j) = fac*sum
    end do

    ! Fold the -0.5*c(1) correction into the coefficients

    cb%c(1) = 0.5_WP*cb%c(1)

    ! Set other components

    cb%x_a = x_a
    cb%x_b = x_b

    cb%n = n

    ! Finish

    return

  end function cheby_t_func_

!****

  $if ($HDF5)

  subroutine read_ (hg, cb)

    type(hgroup_t), intent(inout) :: hg
    type(cheby_t), intent(out)    :: cb

    real(WP)              :: x_a
    real(WP)              :: x_b
    real(WP), allocatable :: c(:)

    ! Read the cheby_t

    call read_attr(hg, 'x_a', x_a)
    call read_attr(hg, 'x_b', x_b)

    call read_dset_alloc(hg, 'c', c)

    cb = cheby_t(x_a, x_b, c)

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

    call write_dset(hg, 'c', cb%c)

    ! Finish

    return

  end subroutine write_

  $endif

!****

  subroutine truncate_ (this, tol)

    class(cheby_t), intent(inout) :: this
    real(WP), intent(in)          :: tol
    
    real(WP) :: c_max
    integer  :: j

    ! Truncate the Chebyshev series to relative tolerance tol (i.e.,
    ! assuming that the discarded terms are all smaller than tol times
    ! the maximum term)

    c_max = MAXVAL(ABS(this%c))

    do j = 1, this%n
       if (ALL(ABS(this%c(j:)) < tol*c_max)) exit
    end do

    if (j <= this%n) then
       this%n = j - 1
       call reallocate(this%c, [this%n])
    end if

    ! Finish

    return

  end subroutine truncate_

!****

  function eval_r_ (this, x) result (f)

    class(cheby_t), intent(in) :: this
    real(WP), intent(in)      :: x
    real(WP)                  :: f

    real(WP) :: d
    real(WP) :: dd
    real(WP) :: y
    real(WP) :: y2
    integer  :: j
    real(WP) :: sv

    ! Evaluate the polynomial at x

    d = 0._WP
    dd = 0._WP

    y = (2._WP*x - this%x_a - this%x_b)/(this%x_b - this%x_a)
    y2 = 2._WP*y

    do j = this%n, 2, -1
       sv = d
       d = y2*d - dd + this%c(j)
       dd = sv
    end do

    f = y*d - dd + this%c(1)

    ! Finish

    return

  end function eval_r_

!****

  function eval_c_ (this, z) result (f)

    class(cheby_t), intent(in) :: this
    complex(WP), intent(in)    :: z
    complex(WP)                :: f

    complex(WP) :: d
    complex(WP) :: dd
    complex(WP) :: y
    complex(WP) :: y2
    integer     :: j
    complex(WP) :: sv

    ! Evaluate the polynomial at x

    d = 0._WP
    dd = 0._WP

    y = (2._WP*z - this%x_a - this%x_b)/(this%x_b - this%x_a)
    y2 = 2._WP*y

    do j = this%n, 2, -1
       sv = d
       d = y2*d - dd + this%c(j)
       dd = sv
    end do

    f = y*d - dd + this%c(1)

    ! Finish

    return

  end function eval_c_

end module gyre_cheby
