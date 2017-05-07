! Module   : gyre_r_ext
! Purpose  : extented-range arithmetic (real)
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
$include 'core_parallel.inc'

module gyre_r_ext

  ! Uses

  use core_kinds
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  real(WP), parameter :: RADIX_WP = REAL(RADIX(1._WP), WP)
  
  ! Derived-type definitions

  $define $OP_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: op_${INFIX}_rx_rx_
    procedure             :: op_${INFIX}_rx_r_
    procedure, pass(that) :: op_${INFIX}_r_rx_
    generic, public       :: operator($OP) => op_${INFIX}_rx_rx_, op_${INFIX}_rx_r_, op_${INFIX}_r_rx_
  $endsub

  type r_ext_t
     private
     real(WP) :: f ! Fractional part
     integer  :: e ! Exponent
   contains
     private
     $OP_DECL(plus,+)
     procedure       :: op_minus_rx_
     generic, public :: operator(-) => op_minus_rx_
     $OP_DECL(minus,-)
     $OP_DECL(times,*)
     $OP_DECL(divide,/)
     $OP_DECL(eq,==)
     $OP_DECL(neq,/=)
     $OP_DECL(lt,<)
     $OP_DECL(gt,>)
     $OP_DECL(le,<=)
     $OP_DECL(ge,>=)
  end type r_ext_t

  ! Interface blocks

  interface r_ext_t
     module procedure r_ext_t_r_
     module procedure r_ext_t_c_
  end interface r_ext_t

  interface real
     module procedure real_
  end interface real

  interface valid
     module procedure valid_
  end interface valid

  interface product
     module procedure product_
  end interface product

  interface exp
     module procedure exp_
  end interface exp

  interface sqrt
     module procedure sqrt_
  end interface sqrt

  interface abs
     module procedure abs_
  end interface abs

  interface fraction
     module procedure fraction_
  end interface fraction

  interface exponent
     module procedure exponent_
  end interface exponent

  interface max
     module procedure max_
  end interface max

  interface min
     module procedure min_
  end interface min

  interface scale
     module procedure scale_
  end interface scale

  $if ($MPI)

  interface send
     module procedure send_0_
     module procedure send_1_
     module procedure send_2_
     module procedure send_3_
     module procedure send_4_
  end interface send

  interface recv
     module procedure recv_0_
     module procedure recv_1_
     module procedure recv_2_
     module procedure recv_3_
     module procedure recv_4_
  end interface recv

  interface recv_any
     module procedure recv_any_0_
     module procedure recv_any_1_
     module procedure recv_any_2_
     module procedure recv_any_3_
     module procedure recv_any_4_
  end interface recv_any

  interface bcast
     module procedure bcast_0_
     module procedure bcast_1_
     module procedure bcast_2_
     module procedure bcast_3_
     module procedure bcast_4_
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_0_
     module procedure bcast_alloc_1_
     module procedure bcast_alloc_2_
     module procedure bcast_alloc_3_
     module procedure bcast_alloc_4_
  end interface bcast_alloc

  interface gatherv
     module procedure gatherv_0_
     module procedure gatherv_1_
     module procedure gatherv_2_
     module procedure gatherv_3_
     module procedure gatherv_4_
  end interface gatherv

  interface allgatherv
     module procedure allgatherv_0_
     module procedure allgatherv_1_
     module procedure allgatherv_2_
     module procedure allgatherv_3_
     module procedure allgatherv_4_
  end interface allgatherv

  $endif

  ! Access specifiers

  private

  public :: r_ext_t
  public :: real
  public :: valid
  public :: product
  public :: abs
  public :: exp
  public :: sqrt
  public :: fraction
  public :: exponent
  public :: max
  public :: min
  public :: scale
  $if ($MPI)
  public :: send
  public :: recv
  public :: recv_any
  public :: bcast
  public :: bcast_alloc
  public :: gatherv
  public :: allgatherv
  $endif

  ! Procedures

contains

  elemental function r_ext_t_r_ (r) result (rx)

    real(WP), intent(in) :: r
    type(r_ext_t)        :: rx

    ! Construct the r_ext_t from the real r

    call split_(r, rx%f, rx%e)

    ! Finish

    return

  end function r_ext_t_r_

  !****

  elemental function r_ext_t_c_ (c) result (rx)

    complex(WP), intent(in) :: c
    type(r_ext_t)           :: rx

    ! Construct the r_ext_t from the complex c

    call split_(REAL(c), rx%f, rx%e)

    ! Finish

    return

  end function r_ext_t_c_

  !****

  elemental function op_plus_rx_rx_ (this, that) result (rx)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(r_ext_t)              :: rx

    real(WP) :: f
    integer  :: e

    ! Evaluate the plus operator

    if (this%f == 0._WP) then

       rx%f = that%f
       rx%e = that%e

    elseif (that%f == 0._WP) then

       rx%f = this%f
       rx%e = this%e

    else

       if (this%e > that%e) then
          f = this%f + real(r_ext_t(that%f, that%e - this%e))
          e = this%e
       else
          f = real(r_ext_t(this%f, this%e - that%e)) + that%f
          e = that%e
       endif

       call split_(f, rx%f, rx%e)
       rx = scale(rx, e)

    endif

    ! Finish

    return

  end function op_plus_rx_rx_

  !****

  elemental function op_minus_rx_ (this) result (rx)

    class(r_ext_t), intent(in) :: this
    type(r_ext_t)              :: rx

    ! Evaluate the unary minus operator

    rx%f = -this%f
    rx%e = this%e

    ! Finish

    return

  end function op_minus_rx_

  !****

  elemental function op_minus_rx_rx_ (this, that) result (rx)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(r_ext_t)              :: rx

    real(WP) :: f
    integer  :: e

    ! Evaluate the minus operator

    if (this%f == 0._WP) then
       
       rx%f = -that%f
       rx%e = that%e

    elseif (that%f == 0._WP) then

       rx%f = this%f
       rx%e = this%e

    else

       if (this%e > that%e) then
          f = this%f - real(r_ext_t(that%f, that%e - this%e))
          e = this%e
       else
          f = real(r_ext_t(this%f, this%e - that%e)) - that%f
          e = that%e
       endif
       
       call split_(f, rx%f, rx%e)
       rx = scale(rx, e)

    endif

    ! Finish

    return

  end function op_minus_rx_rx_

  !****

  elemental function op_times_rx_rx_ (this, that) result (rx)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(r_ext_t)              :: rx

    real(WP) :: f
    integer  :: e

    ! Evaluate the times operator

    if (this%f == 0._WP .OR. that%f == 0._WP) then

       rx = r_ext_t(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_(f, rx%f, rx%e)
       rx = scale(rx, e)

    endif

    ! Finish

    return

  end function op_times_rx_rx_

  !****

  elemental function op_divide_rx_rx_ (this, that) result (rx)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(r_ext_t)              :: rx

    real(WP) :: f
    integer  :: e

    ! Evaluate the divide operator

    if (this%f == 0._WP .AND. that%f /= 0._WP) then

       rx = r_ext_t(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_(f, rx%f, rx%e)
       rx = scale(rx, e)

    endif

    ! Finish

    return

  end function op_divide_rx_rx_

  !****

  elemental function op_eq_rx_rx_ (this, that) result (eq)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: eq

    ! Evaluate the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function op_eq_rx_rx_

  !****

  elemental function op_neq_rx_rx_ (this, that) result (neq)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: neq

    ! Evaluate the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function op_neq_rx_rx_

  !****

  elemental function op_lt_rx_rx_ (this, that) result (lt)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: lt

    real(WP) :: this_s
    real(WP) :: that_s

    ! Evaluate the less-than operator

    if (this%f == 0._WP .OR. that%f == 0._WP) then

       lt = this%f < that%f

    else

       this_s = SIGN(1._WP, this%f)
       that_s = SIGN(1._WP, that%f)

       if (this_s == that_s) then
          if (this_s > 0._WP) then
             lt = (this%e < that%e .OR. (this%e == that%e .AND. this%f < that%f))
          else
             lt = (this%e > that%e .OR. (this%e == that%e .AND. this%f < that%f))
          endif
       else
          lt = this_s < 0._WP
       endif

    endif

    ! Finish

    return

  end function op_lt_rx_rx_

  !****

  elemental function op_gt_rx_rx_ (this, that) result (gt)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: gt

    ! Evaluate the greater-than operator

    gt = that < this

    ! Finish

    return

  end function op_gt_rx_rx_

  !****

  elemental function op_le_rx_rx_ (this, that) result (le)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: le

    ! Evaluate the less-than-or-equal operator

    le = .NOT. this > that

    ! Finish

    return

  end function op_le_rx_rx_
  
  !****

  elemental function op_ge_rx_rx_ (this, that) result (ge)

    class(r_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                    :: ge

    ! Evaluate the greater-than-or-equal operator

    ge = .NOT. this < that

    ! Finish

    return

  end function op_ge_rx_rx_

  !****

  $define $MIXED_OP $sub

  $local $INFIX $1
  $local $OP $2
  $local $TYPE $3

  elemental function op_${INFIX}_rx_r_ (this, that) result ($INFIX)

    class(r_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    $TYPE                      :: $INFIX

    ! Evaluate the operator to mixed ext_real_t/real types

    $INFIX = this $OP r_ext_t(that)

    ! Finish

    return

  end function op_${INFIX}_rx_r_

  elemental function op_${INFIX}_r_rx_ (this, that) result ($INFIX)

    real(WP), intent(in)       :: this
    class(r_ext_t), intent(in) :: that
    $TYPE                      :: $INFIX

    ! Evaluate the operator to mixed real/ext_real_t types

    $INFIX = r_ext_t(this) $OP that

    ! Finish

    return

  end function op_${INFIX}_r_rx_

  $endsub

  $MIXED_OP(plus,+,type(r_ext_t))
  $MIXED_OP(minus,-,type(r_ext_t))
  $MIXED_OP(times,*,type(r_ext_t))
  $MIXED_OP(divide,/,type(r_ext_t))
  $MIXED_OP(lt,<,logical)
  $MIXED_OP(gt,>,logical)
  $MIXED_OP(le,<=,logical)
  $MIXED_OP(ge,>=,logical)
  $MIXED_OP(eq,==,logical)
  $MIXED_OP(neq,/=,logical)

  !****

  elemental function real_ (rx) result (r)

    type(r_ext_t), intent(in) :: rx
    real(WP)                  :: r

    integer :: e_min

    ! Convert the r_ext_t to real

    if (rx%f /= 0._WP) then

       e_min = MINEXPONENT(0._WP)

       if (rx%e >= e_min) then
          r = rx%f*RADIX_WP**rx%e
       else
          r = (rx%f*RADIX_WP**MAX(rx%e-e_min, -DIGITS(0._WP)-1))*RADIX_WP**e_min
       endif

    else

       r = 0._WP

    endif

    ! Finish

    return

  end function real_

  !****

  elemental function valid_ (rx) result (valid)

    type(r_ext_t), intent(in) :: rx
    logical                   :: valid

    ! Determine if rx is valid

    valid = ABS(rx%f) >= RADIX_WP**(-1) .AND. ABS(rx%f) < 1._WP

    ! Finish

    return

  end function valid_

  !****

  function product_ (rx) result (prod_rx)

    type(r_ext_t), intent(in) :: rx(:)
    type(r_ext_t)             :: prod_rx

    integer  :: i
    real(WP) :: g

    ! Calculate the product of the elements of rx

    prod_rx%f = 1._WP
    prod_rx%e = SUM(rx%e)

    do i = 1,SIZE(rx)

       g = prod_rx%f*rx(i)%f
       
       prod_rx%f = FRACTION(g)
       prod_rx = scale(prod_rx, EXPONENT(g))

    end do

    ! Finish

    return

  end function product_

  !****

  elemental function abs_ (rx) result (abs_rx)

    type(r_ext_t), intent(in) :: rx
    type(r_ext_t)             :: abs_rx

    ! Calculate the absolute value of rx

    abs_rx%f = ABS(rx%f)
    abs_rx%e = rx%e

    ! Finish

    return

  end function abs_

  !****

  elemental function exp_ (rx) result (exp_rx)

    type(r_ext_t), intent(in) :: rx
    type(r_ext_t)             :: exp_rx

    real(WP) :: g
    integer  :: e

    ! Calculate the exponential of rx

    g = real(rx)/LOG(RADIX_WP)
    e = FLOOR(g)

    call split_(RADIX_WP**(g-e), exp_rx%f, exp_rx%e)
    exp_rx = scale(exp_rx, e)

    ! Finish

    return

  end function exp_

  !****

  elemental function sqrt_ (rx) result (sqrt_rx)

    type(r_ext_t), intent(in) :: rx
    type(r_ext_t)             :: sqrt_rx

    ! Calculate the square root of rx

    sqrt_rx = r_ext_t(SQRT(rx%f))

    if (MOD(rx%e, 2) == 0) then
       sqrt_rx = scale(sqrt_rx, rx%e/2)
    else
       sqrt_rx = scale(sqrt_rx, (rx%e-1)/2)*SQRT(2._WP)
    endif

    ! Finish

    return

  end function sqrt_

  !****

  elemental function fraction_ (rx) result (fraction_rx)

    type(r_ext_t), intent(in) :: rx
    real(WP)                  :: fraction_rx

    ! Return the fraction part of rx

    fraction_rx = rx%f

    ! Finish

    return

  end function fraction_

  !****

  elemental function exponent_ (rx) result (exponent_rx)

    type(r_ext_t), intent(in) :: rx
    integer                   :: exponent_rx

    ! Return the exponent part of rx

    exponent_rx = rx%e

    ! Finish

    return

  end function exponent_

  !****

  elemental function max_ (rx_a, rx_b) result (max_rx)

    type(r_ext_t), intent(in) :: rx_a
    type(r_ext_t), intent(in) :: rx_b
    type(r_ext_t)             :: max_rx

    real(WP) :: rx_a_s
    real(WP) :: rx_b_s

    ! Return the maximum of rx_a and rx_b

    rx_a_s = SIGN(1._WP, rx_a%f)
    rx_b_s = SIGN(1._WP, rx_b%f)

    if (rx_a_s == rx_b_s) then
       if (rx_a%e > rx_b%e .EQV. rx_a_s > 0._WP) then
          max_rx = rx_a
       elseif (rx_b%e > rx_a%e .EQV. rx_a_s > 0._WP) then
          max_rx = rx_b
       else
          max_rx%f = MAX(rx_a%f, rx_b%f)
          max_rx%e = rx_a%e
       endif
    else
       if (rx_a_s > 0._WP) then
          max_rx = rx_a
       else
          max_rx = rx_b
       endif
    endif

    ! Finish

    return

  end function max_

  !****

  elemental function min_ (rx_a, rx_b) result (min_rx)

    type(r_ext_t), intent(in) :: rx_a
    type(r_ext_t), intent(in) :: rx_b
    type(r_ext_t)             :: min_rx

    ! Return the minimum of rx_a and rx_b

    min_rx = -MAX(-rx_a, -rx_b)

    ! Finish

    return

  end function min_

  !****

  elemental function scale_ (rx, de) result (scale_rx)

    class(r_ext_t), intent(in) :: rx
    integer, intent(in)        :: de
    type(r_ext_t)              :: scale_rx

    ! Scale rx by RADIX_WP**de

    scale_rx%f = rx%f

    if (scale_rx%f /= 0._WP) then
       scale_rx%e = rx%e + de
    else
       scale_rx%e = 0
    endif

    ! Finish

    return

  end function scale_

  !****

  elemental subroutine split_ (x, f, e)

    real(WP), intent(in)  :: x
    real(WP), intent(out) :: f
    integer, intent(out)  :: e

    ! Spit x into fraction and exponent parts

    f = FRACTION(x)
    e = EXPONENT(x)

    ! Finish

    return

  end subroutine split_

  !****

  $if ($MPI)

  $define $SEND $sub

  $local $RX_RANK $1

  subroutine send_${RX_RANK}_ (rx, dest_rank, tag, sync)

    type(r_ext_t), intent(in)     :: rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)           :: dest_rank
    integer, optional, intent(in) :: tag
    logical, optional, intent(in) :: sync

    ! Send the r_ext_t

    call send(rx%f, dest_rank, tag, sync)
    call send(rx%e, dest_rank, tag, sync)

    ! Finish

    return

  end subroutine send_${RX_RANK}_

  $endsub

  $SEND(0)
  $SEND(1)
  $SEND(2)
  $SEND(3)
  $SEND(4)

  !****
  
  $define $RECV $sub

  $local $RX_RANK $1

  subroutine recv_${RX_RANK}_ (rx, src_rank, tag)

    type(r_ext_t), intent(out)    :: rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)           :: src_rank
    integer, optional, intent(in) :: tag

    ! Receive the ext_r_t
    
    call recv(rx%f, src_rank, tag)
    call recv(rx%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_${RX_RANK}_

  $endsub

  $RECV(0)
  $RECV(1)
  $RECV(2)
  $RECV(3)
  $RECV(4)

  !****
  
  $define $RECV_ANY $sub

  $local $RX_RANK $1

  subroutine recv_any_${RX_RANK}_ (rx, src_rank, tag)

    type(r_ext_t), intent(out)    :: rx$ARRAY_SPEC($RX_RANK)
    integer, intent(out)          :: src_rank
    integer, optional, intent(in) :: tag

    ! Receive the r_ext_t

    call recv(rx%f, src_rank, tag)
    call recv(rx%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_any_${RX_RANK}_

  $endsub

  $RECV_ANY(0)
  $RECV_ANY(1)
  $RECV_ANY(2)
  $RECV_ANY(3)
  $RECV_ANY(4)

  !****

  $define $BCAST $sub

  $local $RX_RANK $1

  subroutine bcast_${RX_RANK}_ (rx, root_rank)

    type(r_ext_t), intent(inout) :: rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)          :: root_rank

    ! Broadcast the r_ext_t

    call bcast(rx%f, root_rank)
    call bcast(rx%e, root_rank)

    ! Finish

    return

  end subroutine bcast_${RX_RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)
  $BCAST(2)
  $BCAST(3)
  $BCAST(4)

  !****

  $BCAST_ALLOC(type(r_ext_t),0)
  $BCAST_ALLOC(type(r_ext_t),1)
  $BCAST_ALLOC(type(r_ext_t),2)
  $BCAST_ALLOC(type(r_ext_t),3)
  $BCAST_ALLOC(type(r_ext_t),4)

  !****

  $define $GATHERV $sub

  $local $RX_RANK $1

  subroutine gatherv_${RX_RANK}_ (send_rx, sendcount, recv_rx, recvcounts, displs, root_rank)

    type(r_ext_t), intent(in)    :: send_rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)          :: sendcount
    type(r_ext_t), intent(inout) :: recv_rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)          :: recvcounts(:)
    integer, intent(in)          :: displs(:)
    integer, intent(in)          :: root_rank

    ! Gather the r_ext_t's

    call gatherv(send_rx%f, sendcount, recv_rx%f, recvcounts, displs, root_rank)
    call gatherv(send_rx%e, sendcount, recv_rx%e, recvcounts, displs, root_rank)

    ! Finish

    return

  end subroutine gatherv_${RX_RANK}_

  $endsub

  $GATHERV(0)
  $GATHERV(1)
  $GATHERV(2)
  $GATHERV(3)
  $GATHERV(4)

  !****

  $define $ALLGATHERV $sub

  $local $RX_RANK $1

  subroutine allgatherv_${RX_RANK}_ (rx, recvcounts, displs)

    type(r_ext_t), intent(inout) :: rx$ARRAY_SPEC($RX_RANK)
    integer, intent(in)          :: recvcounts(:)
    integer, intent(in)          :: displs(:)

    ! Gather and share the r_ext_t's

    call allgatherv(rx%f, recvcounts, displs)
    call allgatherv(rx%e, recvcounts, displs)

    ! Finish

    return

  end subroutine allgatherv_${RX_RANK}_

  $endsub

  $ALLGATHERV(0)
  $ALLGATHERV(1)
  $ALLGATHERV(2)
  $ALLGATHERV(3)
  $ALLGATHERV(4)

  $endif

end module gyre_r_ext
