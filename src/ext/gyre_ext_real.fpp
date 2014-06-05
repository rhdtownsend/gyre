! Module   : gyre_ext_real
! Purpose  : extented-range arithmetic (real)
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
$include 'core_parallel.inc'

module gyre_ext_real

  ! Uses

  use core_kinds
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  real(WP), parameter :: RADIX_WP = REAL(RADIX(1._WP), WP)
  
  ! Derived-type definitions

  $define $ARITH_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: arith_${INFIX}_er_er_
    procedure             :: arith_${INFIX}_er_r_
    procedure, pass(that) :: arith_${INFIX}_r_er_
    generic, public       :: operator($OP) => arith_${INFIX}_er_er_, arith_${INFIX}_er_r_, arith_${INFIX}_r_er_
  $endsub

  $define $CMP_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: cmp_${INFIX}_er_er_
    procedure             :: cmp_${INFIX}_er_r_
    procedure, pass(that) :: cmp_${INFIX}_r_er_
    generic, public       :: operator($OP) => cmp_${INFIX}_er_er_, cmp_${INFIX}_er_r_, cmp_${INFIX}_r_er_
  $endsub

  type ext_real_t
     private
     real(WP) :: f ! Fractional part
     integer  :: e ! Exponent
   contains
     private
     $ARITH_DECL(plus,+)
     procedure       :: arith_minus_er_
     generic, public :: operator(-) => arith_minus_er_
     $ARITH_DECL(minus,-)
     $ARITH_DECL(times,*)
     $ARITH_DECL(divide,/)
     $CMP_DECL(eq,==)
     $CMP_DECL(neq,/=)
     $CMP_DECL(lt,<)
     $CMP_DECL(gt,>)
     $CMP_DECL(le,<=)
     $CMP_DECL(ge,>=)
  end type ext_real_t

  ! Interface blocks

  interface ext_real_t
     module procedure ext_real_t_r_
     module procedure ext_real_t_c_
  end interface ext_real_t

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

  public :: ext_real_t
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

  elemental function ext_real_t_r_ (x) result (ex)

    real(WP), intent(in) :: x
    type(ext_real_t)     :: ex

    ! Construct the ext_real_t from the real x

    call split_(x, ex%f, ex%e)

    ! Finish

    return

  end function ext_real_t_r_

!****

  elemental function ext_real_t_c_ (z) result (ex)

    complex(WP), intent(in) :: z
    type(ext_real_t)        :: ex

    ! Construct the ext_real_t from the complex z

    call split_(REAL(z), ex%f, ex%e)

    ! Finish

    return

  end function ext_real_t_c_

!****

  elemental function arith_plus_er_er_ (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the plus operator

    if (this%f == 0._WP) then

       ex%f = that%f
       ex%e = that%e

    elseif (that%f == 0._WP) then

       ex%f = this%f
       ex%e = this%e

    else

       if (this%e > that%e) then
          f = this%f + real(ext_real_t(that%f, that%e - this%e))
          e = this%e
       else
          f = real(ext_real_t(this%f, this%e - that%e)) + that%f
          e = that%e
       endif

       call split_(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_plus_er_er_

!****

  elemental function arith_minus_er_ (this) result (ex)

    class(ext_real_t), intent(in) :: this
    type(ext_real_t)              :: ex

    ! Apply the unary minus operator

    ex%f = -this%f
    ex%e = this%e

    ! Finish

    return

  end function arith_minus_er_

!****

  elemental function arith_minus_er_er_ (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the minus operator

    if (this%f == 0._WP) then
       
       ex%f = -that%f
       ex%e = that%e

    elseif (that%f == 0._WP) then

       ex%f = this%f
       ex%e = this%e

    else

       if (this%e > that%e) then
          f = this%f - real(ext_real_t(that%f, that%e - this%e))
          e = this%e
       else
          f = real(ext_real_t(this%f, this%e - that%e)) - that%f
          e = that%e
       endif
       
       call split_(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_minus_er_er_

!****

  elemental function arith_times_er_er_ (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the times operator

    if (this%f == 0._WP .OR. that%f == 0._WP) then

       ex = ext_real_t(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_times_er_er_

!****

  elemental function arith_divide_er_er_ (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the divide operator

    if (this%f == 0._WP .AND. that%f /= 0._WP) then

       ex = ext_real_t(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_divide_er_er_

!****

  $define $ARITH $sub

  $local $INFIX $1
  $local $OP $2

  elemental function arith_${INFIX}_er_r_ (this, that) result ($INFIX)

    class(ext_real_t), intent(in) :: this
    real(WP), intent(in)          :: that
    type(ext_real_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_real_t(that)

    ! Finish

    return

  end function arith_${INFIX}_er_r_

  elemental function arith_${INFIX}_r_er_ (this, that) result ($INFIX)

    real(WP), intent(in)          :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_real_t(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_r_er_

  $endsub

  $ARITH(plus,+)
  $ARITH(minus,-)
  $ARITH(times,*)
  $ARITH(divide,/)

!****

  elemental function cmp_eq_er_er_ (this, that) result (eq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: eq

    ! Apply the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function cmp_eq_er_er_

!****

  elemental function cmp_neq_er_er_ (this, that) result (neq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: neq

    ! Apply the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function cmp_neq_er_er_

!****

  elemental function cmp_lt_er_er_ (this, that) result (lt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: lt

    real(WP) :: this_s
    real(WP) :: that_s

    ! Apply the less-than operator

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

  end function cmp_lt_er_er_

!****

  elemental function cmp_gt_er_er_ (this, that) result (gt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: gt

    ! Apply the greater-than operator

    gt = that < this

    ! Finish

    return

  end function cmp_gt_er_er_

!****

  elemental function cmp_le_er_er_ (this, that) result (le)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: le

    ! Apply the less-than-or-equal operator

    le = .NOT. this > that

    ! Finish

    return

  end function cmp_le_er_er_
  
!****

  elemental function cmp_ge_er_er_ (this, that) result (ge)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: ge

    ! Apply the greater-than-or-equal operator

    ge = .NOT. this < that

    ! Finish

    return

  end function cmp_ge_er_er_

!****

  $define $CMP $sub

  $local $INFIX $1
  $local $OP $2

  elemental function cmp_${INFIX}_er_r_ (this, that) result ($INFIX)

    class(ext_real_t), intent(in) :: this
    real(WP), intent(in)          :: that
    logical                       :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_real_t(that)

    ! Finish

    return

  end function cmp_${INFIX}_er_r_

  elemental function cmp_${INFIX}_r_er_ (this, that) result ($INFIX)

    real(WP), intent(in)          :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_real_t(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_r_er_

  $endsub

  $CMP(lt,<)
  $CMP(gt,>)
  $CMP(le,<=)
  $CMP(ge,>=)
  $CMP(eq,==)
  $CMP(neq,/=)

!****

  elemental function real_ (ex) result (x)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: x

    integer :: e_min

    ! Convert the ext_real_t to real

    if (ex%f /= 0._WP) then

       e_min = MINEXPONENT(0._WP)

       if (ex%e >= e_min) then
          x = ex%f*RADIX_WP**ex%e
       else
          x = (ex%f*RADIX_WP**MAX(ex%e-e_min, -DIGITS(0._WP)-1))*RADIX_WP**e_min
       endif

    else

       x = 0._WP

    endif

    ! Finish

    return

  end function real_

!****

  elemental function valid_ (ex) result (valid)

    type(ext_real_t), intent(in) :: ex
    logical                      :: valid

    ! Determine if ex is valid

    valid = ABS(ex%f) >= RADIX_WP**(-1) .AND. ABS(ex%f) < 1._WP

    ! Finish

    return

  end function valid_

!****

  function product_ (ex) result (prod_ex)

    type(ext_real_t), intent(in) :: ex(:)
    type(ext_real_t)             :: prod_ex

    integer  :: i
    real(WP) :: g

    ! Calculate the product of the elements of ex

    prod_ex%f = 1._WP
    prod_ex%e = SUM(ex%e)

    do i = 1,SIZE(ex)

       g = prod_ex%f*ex(i)%f
       
       prod_ex%f = FRACTION(g)
       prod_ex = scale(prod_ex, EXPONENT(g))

    end do

    ! Finish

    return

  end function product_

!****

  elemental function abs_ (ex) result (abs_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: abs_ex

    ! Calculate the absolute value of ex

    abs_ex%f = ABS(ex%f)
    abs_ex%e = ex%e

    ! Finish

    return

  end function abs_

!****

  elemental function exp_ (ex) result (exp_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: exp_ex

    real(WP) :: g
    integer  :: e

    ! Calculate the exponential of ex

    g = real(ex)/LOG(RADIX_WP)
    e = FLOOR(g)

    call split_(RADIX_WP**(g-e), exp_ex%f, exp_ex%e)
    exp_ex = scale(exp_ex, e)

    ! Finish

    return

  end function exp_

!****

  elemental function sqrt_ (ex) result (sqrt_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: sqrt_ex

    ! Calculate the square root of ex

    sqrt_ex = ext_real_t(SQRT(ex%f))

    if (MOD(ex%e, 2) == 0) then
       sqrt_ex = scale(sqrt_ex, ex%e/2)
    else
       sqrt_ex = scale(sqrt_ex, (ex%e-1)/2)*SQRT(2._WP)
    endif

    ! Finish

    return

  end function sqrt_

!****

  elemental function fraction_ (ex) result (fraction_ex)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: fraction_ex

    ! Return the fraction part of ex

    fraction_ex = ex%f

    ! Finish

    return

  end function fraction_

!****

  elemental function exponent_ (ex) result (exponent_ex)

    type(ext_real_t), intent(in) :: ex
    integer                      :: exponent_ex

    ! Return the exponent part of ex

    exponent_ex = ex%e

    ! Finish

    return

  end function exponent_

!****

  elemental function max_ (ex_a, ex_b) result (max_ex)

    type(ext_real_t), intent(in) :: ex_a
    type(ext_real_t), intent(in) :: ex_b
    type(ext_real_t)             :: max_ex

    real(WP) :: ex_a_s
    real(WP) :: ex_b_s

    ! Return the maximum of ex_a and ex_b

    ex_a_s = SIGN(1._WP, ex_a%f)
    ex_b_s = SIGN(1._WP, ex_b%f)

    if (ex_a_s == ex_b_s) then
       if (ex_a%e > ex_b%e .EQV. ex_a_s > 0._WP) then
          max_ex = ex_a
       elseif (ex_b%e > ex_a%e .EQV. ex_a_s > 0._WP) then
          max_ex = ex_b
       else
          max_ex%f = MAX(ex_a%f, ex_b%f)
          max_ex%e = ex_a%e
       endif
    else
       if (ex_a_s > 0._WP) then
          max_ex = ex_a
       else
          max_ex = ex_b
       endif
    endif

    ! Finish

    return

  end function max_

!****

  elemental function min_ (ex_a, ex_b) result (min_ex)

    type(ext_real_t), intent(in) :: ex_a
    type(ext_real_t), intent(in) :: ex_b
    type(ext_real_t)             :: min_ex

    ! Return the minimum of ex_a and ex_b

    min_ex = -MAX(-ex_a, -ex_b)

    ! Finish

    return

  end function min_

!****

  elemental function scale_ (ex, de) result (scale_ex)

    class(ext_real_t), intent(in) :: ex
    integer, intent(in)           :: de
    type(ext_real_t)              :: scale_ex

    ! Scale ex by RADIX_WP**de

    scale_ex%f = ex%f
    scale_ex%e = ex%e + de

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

  $local $BUFFER_RANK $1

  subroutine send_${BUFFER_RANK}_ (buffer, dest_rank, tag, sync)

    type(ext_real_t), intent(in)  :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: dest_rank
    integer, optional, intent(in) :: tag
    logical, optional, intent(in) :: sync

    ! Send the buffer

    call send(buffer%f, dest_rank, tag, sync)
    call send(buffer%e, dest_rank, tag, sync)

    ! Finish

    return

  end subroutine send_${BUFFER_RANK}_

  $endsub

  $SEND(0)
  $SEND(1)
  $SEND(2)
  $SEND(3)
  $SEND(4)

!****
  
  $define $RECV $sub

  $local $BUFFER_RANK $1

  subroutine recv_${BUFFER_RANK}_ (buffer, src_rank, tag)

    type(ext_real_t), intent(out) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: src_rank
    integer, optional, intent(in) :: tag

    ! Receive the buffer

    call recv(buffer%f, src_rank, tag)
    call recv(buffer%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_${BUFFER_RANK}_

  $endsub

  $RECV(0)
  $RECV(1)
  $RECV(2)
  $RECV(3)
  $RECV(4)

!****
  
  $define $RECV_ANY $sub

  $local $BUFFER_RANK $1

  subroutine recv_any_${BUFFER_RANK}_ (buffer, src_rank, tag)

    type(ext_real_t), intent(out) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(out)          :: src_rank
    integer, optional, intent(in) :: tag

    ! Receive the buffer

    call recv(buffer%f, src_rank, tag)
    call recv(buffer%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_any_${BUFFER_RANK}_

  $endsub

  $RECV_ANY(0)
  $RECV_ANY(1)
  $RECV_ANY(2)
  $RECV_ANY(3)
  $RECV_ANY(4)

!****

  $define $BCAST $sub

  $local $BUFFER_RANK $1

  subroutine bcast_${BUFFER_RANK}_ (buffer, root_rank)

    type(ext_real_t), intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)             :: root_rank

    ! Broadcast the buffer

    call bcast(buffer%f, root_rank)
    call bcast(buffer%e, root_rank)

    ! Finish

    return

  end subroutine bcast_${BUFFER_RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)
  $BCAST(2)
  $BCAST(3)
  $BCAST(4)

!****

  $BCAST_ALLOC(type(ext_real_t),0)
  $BCAST_ALLOC(type(ext_real_t),1)
  $BCAST_ALLOC(type(ext_real_t),2)
  $BCAST_ALLOC(type(ext_real_t),3)
  $BCAST_ALLOC(type(ext_real_t),4)

!****

  $define $GATHERV $sub

  $local $BUFFER_RANK $1

  subroutine gatherv_${BUFFER_RANK}_ (send_buffer, sendcount, recv_buffer, recvcounts, displs, root_rank)

    type(ext_real_t), intent(in)    :: send_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)             :: sendcount
    type(ext_real_t), intent(inout) :: recv_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)             :: recvcounts(:)
    integer, intent(in)             :: displs(:)
    integer, intent(in)             :: root_rank

    ! Gather the buffers

    call gatherv(send_buffer%f, sendcount, recv_buffer%f, recvcounts, displs, root_rank)
    call gatherv(send_buffer%e, sendcount, recv_buffer%e, recvcounts, displs, root_rank)

    ! Finish

    return

  end subroutine gatherv_${BUFFER_RANK}_

  $endsub

  $GATHERV(0)
  $GATHERV(1)
  $GATHERV(2)
  $GATHERV(3)
  $GATHERV(4)

!****

  $define $ALLGATHERV $sub

  $local $BUFFER_RANK $1

  subroutine allgatherv_${BUFFER_RANK}_ (buffer, recvcounts, displs)

    type(ext_real_t), intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)             :: recvcounts(:)
    integer, intent(in)             :: displs(:)

    ! Gather and share the buffers

    call allgatherv(buffer%f, recvcounts, displs)
    call allgatherv(buffer%e, recvcounts, displs)

    ! Finish

    return

  end subroutine allgatherv_${BUFFER_RANK}_

  $endsub

  $ALLGATHERV(0)
  $ALLGATHERV(1)
  $ALLGATHERV(2)
  $ALLGATHERV(3)
  $ALLGATHERV(4)

  $endif

end module gyre_ext_real
