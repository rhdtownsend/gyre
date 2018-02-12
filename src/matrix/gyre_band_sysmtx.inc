! Incfile  : gyre_band_sysmtx
! Purpose  : system matrix (banded storage, template)
!
! Copyright 2013-2014 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

$if ($DOUBLE_PRECISION)
$if ($T eq 'r')
$local $X D
$else
$local $X Z
$endif
$else
$if ($T eq 'r')
$local $X S
$else
$local $X C
$endif
$endif

module gyre_${T}_band_sysmtx

  ! Uses

  use core_kinds
  use core_parallel
  use core_linalg

  use gyre_ext
  use gyre_linalg
  use gyre_sysmtx

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: ALGO_LEN = 16

  ! Derived-type definitions

  type, extends (${T}_sysmtx_t) :: ${T}_band_sysmtx_t
     private
     $TYPE(WP), allocatable        :: A_b(:,:) ! Banded matrix
     integer, allocatable          :: ipiv(:)  ! Pivot indices
     $TYPE(WP), allocatable        :: scl_i(:) ! Inner boundary scales
     $TYPE(WP), allocatable        :: scl_o(:) ! Outer boundary scales
     type(${T}_ext_t), allocatable :: scl(:)   ! Block scales
     integer                       :: n_ul     ! Number of sub-/super-diagonals
   contains
     private
     procedure, public :: set_B_i
     procedure, public :: set_B_o
     procedure, public :: set_E
     procedure         :: set_row_
     procedure         :: get_row_
     procedure, public :: factor
     procedure         :: scale_rows_
     procedure, public :: det
     procedure, public :: soln_vec_hom
     procedure, public :: soln_vec_inhom
     procedure, public :: resd_vec
  end type ${T}_band_sysmtx_t

  ! Interfaces

  interface ${T}_band_sysmtx_t
     module procedure ${T}_band_sysmtx_t_
  end interface ${T}_band_sysmtx_t

  ! Access specifiers

  private

  public :: ${T}_band_sysmtx_t

  ! Procedures

contains

  function ${T}_band_sysmtx_t_ (n, n_e, n_i, n_o) result (sm)

    integer, intent(in)      :: n
    integer, intent(in)      :: n_e
    integer, intent(in)      :: n_i
    integer, intent(in)      :: n_o
    type(${T}_band_sysmtx_t) :: sm

    ! Construct the sysmtx_t

    ! Note that an additional n_ul rows are added to A_b to provide
    ! space for fill-in during factorization

    sm%n_ul = n_e + n_i - 1

    allocate(sm%A_b(3*sm%n_ul+1,n_e*(n+1)))
    allocate(sm%ipiv(n_e*(n+1)))

    allocate(sm%scl_i(n_i))
    allocate(sm%scl_o(n_o))

    allocate(sm%scl(n))

    sm%n = n
    sm%n_e = n_e
    sm%n_i = n_i
    sm%n_o = n_o

    ! Finish

    return

  end function ${T}_band_sysmtx_t_

  !****

  subroutine set_B_i (this, B, scl)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                    :: B(:,:)
    $TYPE(WP), intent(in)                    :: scl(:)

    integer :: i_0
    integer :: j_0
    integer :: i
    
    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Set the inner boundary conditions

    i_0 = 1
    j_0 = 1

    do i = 1, this%n_i
       call this%set_row_(i_0+i-1, j_0, B(i,:))
    end do

    this%scl_i = scl

    ! Finish

    return

  end subroutine set_B_i

  !****

  subroutine set_B_o (this, B, scl)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                    :: B(:,:)
    $TYPE(WP), intent(in)                    :: scl(:)

    integer :: i_0
    integer :: j_0
    integer :: i

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Set the outer boundary conditions

    i_0 = this%n*this%n_e + this%n_i + 1
    j_0 = this%n*this%n_e + 1

    do i = 1, this%n_o
       call this%set_row_(i_0+i-1, j_0, B(i,:))
    end do

    this%scl_o = scl

    ! Finish

    return

  end subroutine set_B_o

  !****

  subroutine set_E (this, k, E_l, E_r, scl)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    integer, intent(in)                      :: k
    $TYPE(WP), intent(in)                    :: E_l(:,:)
    $TYPE(WP), intent(in)                    :: E_r(:,:)
    type(${T}_ext_t), intent(in)             :: scl

    integer   :: i_0
    integer   :: j_0
    integer   :: i
    $TYPE(WP) :: R(2*this%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    $ASSERT(k >= 1,Invalid block index)
    $ASSERT(k <= this%n,Invalid block index)

    ! Set the block

    i_0 = this%n_e*(k-1) + this%n_i + 1
    j_0 = this%n_e*(k-1) + 1

    do i = 1, this%n_e

       R = [E_l(i,:),E_r(i,:)]

       call this%set_row_(i_0+i-1, j_0, R)

    end do

    this%scl(k) = scl

    ! Finish

    return

  end subroutine set_E

  !****

  subroutine set_row_ (this, i_0, j_0, R)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    integer, intent(in)                      :: i_0
    integer, intent(in)                      :: j_0
    $TYPE(WP), intent(in)                    :: R(:)

    integer :: j
    integer :: i_b
    integer :: j_b

    $ASSERT_DEBUG(ABS(i_0-j_0) <= this%n_ul,Invalid indices)
    $ASSERT_DEBUG(ABS(i_0-(j_0+SIZE(R)-1)) <= this%n_ul,Invalid indices)

    ! Set data in row i_0, starting at column j_0. The rest of the row
    ! is zeroed out

    do j = MAX(i_0-this%n_ul, 1), j_0-1
       i_b = 2*this%n_ul + 1 + i_0 - j
       j_b = j
       this%A_b(i_b, j_b) = 0._WP
    end do

    do j = j_0, j_0+SIZE(R)-1
       i_b = 2*this%n_ul + 1 + i_0 - j
       j_b = j
       this%A_b(i_b, j_b) = R(j-j_0+1)
    end do

    do j = j_0+SIZE(R), MIN(i_0+this%n_ul, this%n_e*(this%n+1))
       i_b = 2*this%n_ul + 1 + i_0 - j
       j_b = j
       this%A_b(i_b, j_b) = 0._WP
    end do

    ! Finish

    return

  end subroutine set_row_

  !****

  subroutine get_row_ (this, i_0, j_0, R)

    class(${T}_band_sysmtx_t), intent(in) :: this
    integer, intent(in)                   :: i_0
    integer, intent(in)                   :: j_0
    $TYPE(WP), intent(out)                :: R(:)

    integer :: j
    integer :: i_b
    integer :: j_b

    $ASSERT_DEBUG(ABS(i_0-j_0) <= this%n_ul,Invalid indices)
    $ASSERT_DEBUG(ABS(i_0-(j_0+SIZE(R)-1)) <= this%n_ul,Invalid indices)

    ! Get data from row i_0, starting at column j_0

    do j = j_0, j_0+SIZE(R)-1
       i_b = 2*this%n_ul + 1 + i_0 - j
       j_b = j
       R(j-j_0+1) = this%A_b(i_b, j_b)
    end do

    ! Finish

    return

  end subroutine get_row_

  !****

  subroutine factor (this)

    class(${T}_band_sysmtx_t), intent(inout) :: this

    $TYPE(WP), parameter :: ONE = 1._WP

    integer :: n
    integer :: m
    integer :: info

    ! Factorize the sysmtx using LU decomposition

    call this%scale_rows_()

    n = SIZE(this%A_b, 1)
    m = SIZE(this%A_b, 2)

    call ${X}GBTRF(m, m, this%n_ul, this%n_ul, this%A_b, n, this%ipiv, info)
    $ASSERT(info == 0 .OR. info > m-this%n_e,Non-zero return from XGBTRF)

    ! Finish

    return

  end subroutine factor

  !****

  subroutine scale_rows_ (this)

    class(${T}_band_sysmtx_t), intent(inout) :: this

    integer  :: i
    integer  :: j
    real(WP) :: scl
    integer  :: i_b
    integer  :: j_b

    ! Scale the rows of the sysmtx to have maximum absolute value of unity

    do i = 1, this%n_e*(this%n+1)

       scl = 0._WP

       do j = MAX(i-this%n_ul, 1), MIN(i+this%n_ul, this%n_e*(this%n+1))
          i_b = 2*this%n_ul + 1 + i - j
          j_b = j
          scl = MAX(scl, ABS(this%A_b(i_b,j_b)))
       end do

       do j = MAX(i-this%n_ul, 1), MIN(i+this%n_ul, this%n_e*(this%n+1))
          i_b = 2*this%n_ul + 1 + i - j
          j_b = j
          this%A_b(i_b,j_b) = this%A_b(i_b,j_b)/scl
       end do

       if (i <= this%n_i) then

          this%scl_i(i) = this%scl_i(i)*scl

       elseif (i > this%n_e*this%n + this%n_i) then

          associate (i_ => i-this%n_e*this%n-this%n_i)
            this%scl_o(i_) = this%scl_o(i_)*scl
          end associate

       else

          associate (k => (i-1+this%n_o)/this%n_e)
            this%scl(k) = this%scl(k)*scl
          end associate

       endif

    end do

    ! Finish

    return

  end subroutine scale_rows_

  !****

  function det (this)

    class(${T}_band_sysmtx_t), intent(in) :: this
    type(${T}_ext_t)                      :: det

    integer :: j

    ! Evaluate the determinant

    det = product([${T}_ext_t(this%A_b(2*this%n_ul+1,:)),${T}_ext_t(this%scl_i),this%scl,${T}_ext_t(this%scl_o)])
 
    do j = 1, SIZE(this%ipiv)
       if(this%ipiv(j) /= j) det = -det
    enddo

    ! Finish

    return

  end function det

  !****

  function soln_vec_hom (this) result (v)
  
    class(${T}_band_sysmtx_t), intent(in) :: this
    $TYPE(WP)                             :: v(this%n_e*(this%n+1))

    integer                :: i_s
    real(WP)               :: A_s
    integer                :: i
    $TYPE(WP), allocatable :: A_b(:,:)
    $TYPE(WP), allocatable :: B(:,:)
    integer, allocatable   :: ipiv(:)
    integer                :: n2
    integer                :: j
    integer                :: info

    ! Evaluate the solution vector v of the homogeneous linear system
    ! S v = 0. It is assumed that the nullity nul(S) >= 1

    associate (n => this%n, n_e => this%n_e, n_ul => this%n_ul)

      ! Locate the smallest element on the diagonal of the outer
      ! block (this will be taken to be the singular element)

      i_s = 0
      A_s = HUGE(0._WP)

      sing_loop : do i = n_e*n+1, n_e*(n+1)
         if (ABS(this%A_b(2*n_ul+1,i)) < A_s) then
            A_s = ABS(this%A_b(2*n_ul+1,i))
            i_s = i
         end if
      end do sing_loop

      ! Set up the reduced banded system

      allocate(A_b(3*n_ul+1,i_s-1))
      allocate(B(i_s-1,1))
      
      allocate(ipiv(i_s-1))

      A_b(:2*n_ul+1,:) = this%A_b(:2*n_ul+1,:i_s-1)
      A_b(2*n_ul+2:,:) = 0._WP

      n2 = MIN(2*n_ul, i_s-1)

      B(:i_s-n2-1,1) = 0._WP
      B(i_s-n2:,1) = -this%A_b(2*n_ul+1-n2:2*n_ul,i_s)

      do j = 1,i_s-1
         ipiv(j) = j
      enddo

      ! Solve for the 1:i_s-1 components of the null vector

      call ${X}GBTRS('N', SIZE(A_b, 2), n_ul, n_ul, 1, A_b, SIZE(A_b, 1), ipiv, B, SIZE(B, 1), info)
      $ASSERT(info == 0,Non-zero return from XGBTRS)

      v(:i_s-1) = B(:,1)

    end associate
       
    ! Fill in the remaining components of the null vector
    
    v(i_s) = 1._WP
    v(i_s+1:) = 0._WP

    ! Finish

    return

  end function soln_vec_hom

  !****

  function soln_vec_inhom (this, w_i, w_o) result (v)
  
    class(${T}_band_sysmtx_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: w_i(:)
    $TYPE(WP), intent(in)                 :: w_o(:)
    $TYPE(WP)                             :: v(this%n_e*(this%n+1))

    $TYPE(WP) :: B(this%n_e*(this%n+1),1)
    integer   :: info

    $CHECK_BOUNDS(SIZE(w_i),this%n_i)
    $CHECK_BOUNDS(SIZE(w_o),this%n_o)

    ! Evaluate the solution vector v of the inhomogeneous linear
    ! system S v = w. It is assumed that the right-hand side vector w
    ! has non-zero components in only the n_i first and n_o last rows
    ! (corresponding to the inner and outer boundary
    ! conditions). These components are supplied in w_i and w_o,
    ! respectively.
    
    associate (n => this%n, n_e => this%n_e, n_i => this%n_i, n_ul => this%n_ul)

      ! Set up the banded system

      B(:n_i,1) = w_i/this%scl_i
      B(n_i+1:n_i+n*n_e,1) = 0._WP
      B(n_i+n*n_e+1:,1) = w_o/this%scl_o

      ! Solve for v

      call ${X}GBTRS('N', n_e*(n+1), n_ul, n_ul, 1, &
                     this%A_b, SIZE(this%A_b, 1), this%ipiv, B, n_e*(n+1), info)
      $ASSERT(info == 0,Non-zero return from XGBTRS)

      v = B(:,1)

    end associate
       
    ! Finish

    return

  end function soln_vec_inhom

  !****

  function resd_vec (this, v) result (w)

    class(${T}_band_sysmtx_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: v(:)
    $TYPE(WP)                             :: w(this%n_e*(this%n+1))

    integer   :: i
    integer   :: j_a
    integer   :: j_b
    $TYPE(WP) :: R(this%n_e*(this%n+1))

    $CHECK_BOUNDS(SIZE(v),this%n_e*(this%n+1))

    ! Evaluate the residual vector w = S v

    associate (n => this%n, n_e => this%n_e)

      !$OMP PARALLEL DO PRIVATE (j_a, j_b, R)
      multiply_loop : do i = 1, n_e*(n+1)

         j_a = MAX(i-this%n_ul, 1)
         j_b = MIN(i+this%n_ul, n_e*(n+1))

         call this%get_row_(i, j_a, R(j_a:j_b))

         w(i) = SUM(R(j_a:j_b)*v(j_a:j_b))

      end do multiply_loop

    end associate

    ! Finish

    return

  end function resd_vec
  
end module gyre_${T}_band_sysmtx

$endsub

