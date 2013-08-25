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
    procedure             :: arith_${INFIX}_er_er
    procedure             :: arith_${INFIX}_er_r
    procedure, pass(that) :: arith_${INFIX}_r_er
    generic, public       :: operator($OP) => arith_${INFIX}_er_er, arith_${INFIX}_er_r, arith_${INFIX}_r_er
  $endsub

  $define $CMP_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: cmp_${INFIX}_er_er
    procedure             :: cmp_${INFIX}_er_r
    procedure, pass(that) :: cmp_${INFIX}_r_er
    generic, public       :: operator($OP) => cmp_${INFIX}_er_er, cmp_${INFIX}_er_r, cmp_${INFIX}_r_er
  $endsub

  type ext_real_t
     private
     real(WP) :: f ! Fractional part
     integer  :: e ! Exponent
   contains
     private
     $ARITH_DECL(plus,+)
     procedure       :: arith_minus_er
     generic, public :: operator(-) => arith_minus_er
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

  $if($MPI)

  interface send
     module procedure send_er_0
     module procedure send_er_1
     module procedure send_er_2
     module procedure send_er_3
     module procedure send_er_4
  end interface send

  interface recv
     module procedure recv_er_0
     module procedure recv_er_1
     module procedure recv_er_2
     module procedure recv_er_3
     module procedure recv_er_4
  end interface recv

  interface recv_any
     module procedure recv_any_er_0
     module procedure recv_any_er_1
     module procedure recv_any_er_2
     module procedure recv_any_er_3
     module procedure recv_any_er_4
  end interface recv_any

  interface bcast
     module procedure bcast_er_0
     module procedure bcast_er_1
     module procedure bcast_er_2
     module procedure bcast_er_3
     module procedure bcast_er_4
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_er_0
     module procedure bcast_alloc_er_1
     module procedure bcast_alloc_er_2
     module procedure bcast_alloc_er_3
     module procedure bcast_alloc_er_4
  end interface bcast_alloc

  interface gatherv
     module procedure gatherv_er_0
     module procedure gatherv_er_1
     module procedure gatherv_er_2
     module procedure gatherv_er_3
     module procedure gatherv_er_4
  end interface gatherv

  interface allgatherv
     module procedure allgatherv_er_0
     module procedure allgatherv_er_1
     module procedure allgatherv_er_2
     module procedure allgatherv_er_3
     module procedure allgatherv_er_4
  end interface allgatherv

  $endif

  interface ext_real
     module procedure ext_real_r
     module procedure ext_real_c
  end interface ext_real

  interface real
     module procedure real_er
  end interface real

  interface valid
     module procedure valid_er
  end interface valid

  interface product
     module procedure product_er
  end interface product

  interface exp
     module procedure exp_er
  end interface exp

  interface abs
     module procedure abs_er
  end interface abs

  interface fraction
     module procedure fraction_er
  end interface fraction

  interface exponent
     module procedure exponent_er
  end interface exponent

  interface max
     module procedure max_er
  end interface max

  interface min
     module procedure min_er
  end interface min

  interface scale
     module procedure scale_er
  end interface scale

  ! Access specifiers

  private

  $if($MPI)
  public :: send
  public :: recv
  public :: recv_any
  public :: bcast
  public :: bcast_alloc
  public :: gatherv
  public :: allgatherv
  $endif
  public :: ext_real_t
  public :: ext_real
  public :: real
  public :: valid
  public :: product
  public :: abs
  public :: exp
  public :: fraction
  public :: exponent
  public :: max
  public :: min
  public :: scale

  ! Procedures

contains

  $if($MPI)

  $define $SEND $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine send_${INFIX}_${BUFFER_RANK} (buffer, dest_rank, tag, sync)

    $BUFFER_TYPE, intent(in)      :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: dest_rank
    integer, intent(in), optional :: tag
    logical, intent(in), optional :: sync

    ! Send the buffer

    call send(buffer%f, dest_rank, tag, sync)
    call send(buffer%e, dest_rank, tag, sync)

    ! Finish

    return

  end subroutine send_${INFIX}_${BUFFER_RANK}

  $endsub

  $SEND(er,type(ext_real_t),0)
  $SEND(er,type(ext_real_t),1)
  $SEND(er,type(ext_real_t),2)
  $SEND(er,type(ext_real_t),3)
  $SEND(er,type(ext_real_t),4)

  $SEND(ec,type(ext_complex_t),0)
  $SEND(ec,type(ext_complex_t),1)
  $SEND(ec,type(ext_complex_t),2)
  $SEND(ec,type(ext_complex_t),3)
  $SEND(ec,type(ext_complex_t),4)

!****
  
  $define $RECV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine recv_${INFIX}_${BUFFER_RANK} (buffer, src_rank, tag)

    $BUFFER_TYPE, intent(out)     :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: src_rank
    integer, intent(in), optional :: tag

    ! Receive the buffer

    call recv(buffer%f, src_rank, tag)
    call recv(buffer%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_${INFIX}_${BUFFER_RANK}

  $endsub

  $RECV(er,type(ext_real_t),0)
  $RECV(er,type(ext_real_t),1)
  $RECV(er,type(ext_real_t),2)
  $RECV(er,type(ext_real_t),3)
  $RECV(er,type(ext_real_t),4)

  $RECV(ec,type(ext_complex_t),0)
  $RECV(ec,type(ext_complex_t),1)
  $RECV(ec,type(ext_complex_t),2)
  $RECV(ec,type(ext_complex_t),3)
  $RECV(ec,type(ext_complex_t),4)
  
!****
  
  $define $RECV_ANY $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine recv_any_${INFIX}_${BUFFER_RANK} (buffer, src_rank, tag)

    $BUFFER_TYPE, intent(out)     :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(out)          :: src_rank
    integer, intent(in), optional :: tag

    ! Receive the buffer

    call recv(buffer%f, src_rank, tag)
    call recv(buffer%e, src_rank, tag)

    ! Finish

    return

  end subroutine recv_any_${INFIX}_${BUFFER_RANK}

  $endsub

  $RECV_ANY(er,type(ext_real_t),0)
  $RECV_ANY(er,type(ext_real_t),1)
  $RECV_ANY(er,type(ext_real_t),2)
  $RECV_ANY(er,type(ext_real_t),3)
  $RECV_ANY(er,type(ext_real_t),4)

  $RECV_ANY(ec,type(ext_complex_t),0)
  $RECV_ANY(ec,type(ext_complex_t),1)
  $RECV_ANY(ec,type(ext_complex_t),2)
  $RECV_ANY(ec,type(ext_complex_t),3)
  $RECV_ANY(ec,type(ext_complex_t),4)
  
!****

  $define $BCAST $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine bcast_${INFIX}_${BUFFER_RANK} (buffer, root_rank)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: root_rank

    ! Broadcast the buffer

    call bcast(buffer%f, root_rank)
    call bcast(buffer%e, root_rank)

    ! Finish

    return

  end subroutine bcast_${INFIX}_${BUFFER_RANK}

  $endsub

  $BCAST(er,type(ext_real_t),0)
  $BCAST(er,type(ext_real_t),1)
  $BCAST(er,type(ext_real_t),2)
  $BCAST(er,type(ext_real_t),3)
  $BCAST(er,type(ext_real_t),4)

  $BCAST(ec,type(ext_complex_t),0)
  $BCAST(ec,type(ext_complex_t),1)
  $BCAST(ec,type(ext_complex_t),2)
  $BCAST(ec,type(ext_complex_t),3)
  $BCAST(ec,type(ext_complex_t),4)
  
!****

  $BCAST_ALLOC(er,type(ext_real_t),0)
  $BCAST_ALLOC(er,type(ext_real_t),1)
  $BCAST_ALLOC(er,type(ext_real_t),2)
  $BCAST_ALLOC(er,type(ext_real_t),3)
  $BCAST_ALLOC(er,type(ext_real_t),4)

  $BCAST_ALLOC(ec,type(ext_complex_t),0)
  $BCAST_ALLOC(ec,type(ext_complex_t),1)
  $BCAST_ALLOC(ec,type(ext_complex_t),2)
  $BCAST_ALLOC(ec,type(ext_complex_t),3)
  $BCAST_ALLOC(ec,type(ext_complex_t),4)

!****

  $define $GATHERV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $TYPE $3
  $local $BUFFER_RANK $4

  subroutine gatherv_${INFIX}_${BUFFER_RANK} (send_buffer, sendcount, recv_buffer, recvcounts, displs, root_rank)

    $BUFFER_TYPE, intent(in)    :: send_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: sendcount
    $BUFFER_TYPE, intent(inout) :: recv_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: recvcounts(:)
    integer, intent(in)         :: displs(:)
    integer, intent(in)         :: root_rank

    ! Gather the buffers

    call gatherv(send_buffer%f, sendcount, recv_buffer%f, recvcounts, displs, root_rank)
    call gatherv(send_buffer%e, sendcount, recv_buffer%e, recvcounts, displs, root_rank)

    ! Finish

    return

  end subroutine gatherv_${INFIX}_${BUFFER_RANK}

  $endsub

  $GATHERV(er,type(ext_real_t),real,0)
  $GATHERV(er,type(ext_real_t),real,1)
  $GATHERV(er,type(ext_real_t),real,2)
  $GATHERV(er,type(ext_real_t),real,3)
  $GATHERV(er,type(ext_real_t),real,4)

  $GATHERV(ec,type(ext_complex_t),complex,0)
  $GATHERV(ec,type(ext_complex_t),complex,1)
  $GATHERV(ec,type(ext_complex_t),complex,2)
  $GATHERV(ec,type(ext_complex_t),complex,3)
  $GATHERV(ec,type(ext_complex_t),complex,4)
  
!****

  $define $ALLGATHERV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine allgatherv_${INFIX}_${BUFFER_RANK} (buffer, recvcounts, displs)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: recvcounts(:)
    integer, intent(in)         :: displs(:)

    ! Gather and share the buffers

    call allgatherv(buffer%f, recvcounts, displs)
    call allgatherv(buffer%e, recvcounts, displs)

    ! Finish

    return

  end subroutine allgatherv_${INFIX}_${BUFFER_RANK}

  $endsub

  $ALLGATHERV(er,type(ext_real_t),0)
  $ALLGATHERV(er,type(ext_real_t),1)
  $ALLGATHERV(er,type(ext_real_t),2)
  $ALLGATHERV(er,type(ext_real_t),3)
  $ALLGATHERV(er,type(ext_real_t),4)

  $ALLGATHERV(ec,type(ext_complex_t),0)
  $ALLGATHERV(ec,type(ext_complex_t),1)
  $ALLGATHERV(ec,type(ext_complex_t),2)
  $ALLGATHERV(ec,type(ext_complex_t),3)
  $ALLGATHERV(ec,type(ext_complex_t),4)
  
  $endif

!****

  elemental function arith_plus_er_er (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the op plus operator

    if(this%f == 0._WP) then

       ex%f = that%f
       ex%e = that%e

    elseif(that%f == 0._WP) then

       ex%f = this%f
       ex%e = this%e

    else

       if(this%e > that%e) then
          f = this%f + that%f*RADIX_WP**(that%e - this%e)
          e = this%e
       else
          f = this%f*RADIX_WP**(this%e - that%e) + that%f
          e = that%e
       endif

       call split_r(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_plus_er_er

!****

  elemental function arith_minus_er (this) result (ex)

    class(ext_real_t), intent(in) :: this
    type(ext_real_t)              :: ex

    ! Apply the unary minus operator

    ex%f = -this%f
    ex%e = this%e

    ! Finish

    return

  end function arith_minus_er

!****

  elemental function arith_minus_er_er (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the minus operator

    if(this%f == 0._WP) then
       
       ex%f = -that%f
       ex%e = that%e

    elseif(that%f == 0._WP) then

       ex%f = this%f
       ex%e = this%e

    else

       if(this%e > that%e) then
          f = this%f - that%f*RADIX_WP**(that%e - this%e)
          e = this%e
       else
          f = this%f*RADIX_WP**(this%e - that%e) - that%f
          e = that%e
       endif
       
       call split_r(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_minus_er_er

!****

  elemental function arith_times_er_er (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the times operator

    if(this%f == 0._WP .OR. that%f == 0._WP) then

       ex = ext_real(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_r(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_times_er_er

!****

  elemental function arith_divide_er_er (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    real(WP) :: f
    integer  :: e

    ! Apply the divide operator

    if(this%f == 0._WP .AND. that%f /= 0._WP) then

       ex = ext_real(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_r(f, ex%f, ex%e)
       ex = scale(ex, e)

    endif

    ! Finish

    return

  end function arith_divide_er_er

!****

  $define $ARITH $sub

  $local $INFIX $1
  $local $OP $2

  elemental function arith_${INFIX}_er_r (this, that) result ($INFIX)

    class(ext_real_t), intent(in) :: this
    real(WP), intent(in)          :: that
    type(ext_real_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_real(that)

    ! Finish

    return

  end function arith_${INFIX}_er_r

  elemental function arith_${INFIX}_r_er (this, that) result ($INFIX)

    real(WP), intent(in)          :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_real(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_r_er

  $endsub

  $ARITH(plus,+)
  $ARITH(minus,-)
  $ARITH(times,*)
  $ARITH(divide,/)

!****

  elemental function cmp_eq_er_er (this, that) result (eq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: eq

    ! Apply the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function cmp_eq_er_er

!****

  elemental function cmp_neq_er_er (this, that) result (neq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: neq

    ! Apply the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function cmp_neq_er_er

!****

  elemental function cmp_lt_er_er (this, that) result (lt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: lt

    real(WP) :: this_s
    real(WP) :: that_s

    ! Apply the less-than operator

    if(this%f == 0._WP .OR. that%f == 0._WP) then

       lt = this%f < that%f

    else

       this_s = SIGN(1._WP, this%f)
       that_s = SIGN(1._WP, that%f)

       if(this_s == that_s) then
          if(this_s > 0._WP) then
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

  end function cmp_lt_er_er

!****

  elemental function cmp_gt_er_er (this, that) result (gt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: gt

    ! Apply the greater-than operator

    gt = that < this

    ! Finish

    return

  end function cmp_gt_er_er

!****

  elemental function cmp_le_er_er (this, that) result (le)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: le

    ! Apply the less-than-or-equal operator

    le = .NOT. this > that

    ! Finish

    return

  end function cmp_le_er_er
  
!****

  elemental function cmp_ge_er_er (this, that) result (ge)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: ge

    ! Apply the greater-than-or-equal operator

    ge = .NOT. this < that

    ! Finish

    return

  end function cmp_ge_er_er

!****

  $define $CMP $sub

  $local $INFIX $1
  $local $OP $2

  elemental function cmp_${INFIX}_er_r (this, that) result ($INFIX)

    class(ext_real_t), intent(in) :: this
    real(WP), intent(in)          :: that
    logical                       :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_real(that)

    ! Finish

    return

  end function cmp_${INFIX}_er_r

  elemental function cmp_${INFIX}_r_er (this, that) result ($INFIX)

    real(WP), intent(in)          :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_real(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_r_er

  $endsub

  $CMP(lt,<)
  $CMP(gt,>)
  $CMP(le,<=)
  $CMP(ge,>=)
  $CMP(eq,==)
  $CMP(neq,/=)

!****

  elemental function ext_real_r (x) result (ex)

    real(WP), intent(in) :: x
    type(ext_real_t)     :: ex

    ! Convert real to ext_real

    call split_r(x, ex%f, ex%e)

    ! Finish

    return

  end function ext_real_r

!****

  elemental function ext_real_c (z) result (ex)

    complex(WP), intent(in) :: z
    type(ext_real_t)        :: ex

    ! Convert complex to ext_real

    call split_r(REAL(z), ex%f, ex%e)

    ! Finish

    return

  end function ext_real_c

!****

  elemental function real_er (ex) result (x)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: x

    ! Convert ext_real to real

    if(ex%f /= 0._WP) then
       x = ex%f*RADIX_WP**ex%e
    else
       x = 0._WP
    endif

    ! Finish

    return

  end function real_er

!****

  elemental function valid_er (ex) result (valid)

    type(ext_real_t), intent(in) :: ex
    logical                      :: valid

    ! Determine if ext_real is valid

    valid = ABS(ex%f) >= RADIX_WP**(-1) .AND. ABS(ex%f) < 1._WP

    ! Finish

    return

  end function valid_er

!****

  function product_er (ex) result (prod_ex)

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

  end function product_er

!****

  elemental function abs_er (ex) result (abs_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: abs_ex

    ! Calculate the absolute value of ex

    abs_ex%f = ABS(ex%f)
    abs_ex%e = ex%e

    ! Finish

    return

  end function abs_er

!****

  elemental function exp_er (ex) result (exp_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: exp_ex

    real(WP) :: g
    integer  :: e

    ! Calculate the exponential of ex

    g = real(ex)/LOG(RADIX_WP)
    e = FLOOR(g)

    call split_r(RADIX_WP**(g-e), exp_ex%f, exp_ex%e)
    exp_ex = scale(exp_ex, e)

    ! Finish

    return

  end function exp_er

!****

  elemental function fraction_er (ex) result (fraction_ex)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: fraction_ex

    ! Return the fraction part of ex

    fraction_ex = ex%f

    ! Finish

    return

  end function fraction_er

!****

  elemental function exponent_er (ex) result (exponent_ex)

    type(ext_real_t), intent(in) :: ex
    integer                      :: exponent_ex

    ! Return the exponent part of ex

    exponent_ex = ex%e

    ! Finish

    return

  end function exponent_er

!****

  elemental function max_er (ex_a, ex_b) result (max_ex)

    type(ext_real_t), intent(in) :: ex_a
    type(ext_real_t), intent(in) :: ex_b
    type(ext_real_t)             :: max_ex

    real(WP) :: ex_a_s
    real(WP) :: ex_b_s

    ! Return the maximum of ex_a and ex_b

    ex_a_s = SIGN(1._WP, ex_a%f)
    ex_b_s = SIGN(1._WP, ex_b%f)

    if(ex_a_s == ex_b_s) then
       if(ex_a%e > ex_b%e .EQV. ex_a_s > 0._WP) then
          max_ex = ex_a
       elseif(ex_b%e > ex_a%e .EQV. ex_a_s > 0._WP) then
          max_ex = ex_b
       else
          max_ex%f = MAX(ex_a%f, ex_b%f)
          max_ex%e = ex_a%e
       endif
    else
       if(ex_a_s > 0._WP) then
          max_ex = ex_a
       else
          max_ex = ex_b
       endif
    endif

    ! Finish

    return

  end function max_er

!****

  elemental function min_er (ex_a, ex_b) result (min_ex)

    type(ext_real_t), intent(in) :: ex_a
    type(ext_real_t), intent(in) :: ex_b
    type(ext_real_t)             :: min_ex

    ! Return the minimum of ex_a and ex_b

    min_ex = -MAX(-ex_a, -ex_b)

    ! Finish

    return

  end function min_er

!****

  elemental function scale_er (ex, de) result (scale_ex)

    class(ext_real_t), intent(in) :: ex
    integer, intent(in)           :: de
    type(ext_real_t)              :: scale_ex

    ! Scale ex by RADIX_WP**de

    scale_ex%f = ex%f
    scale_ex%e = ex%e + de

    ! Finish

    return

  end function scale_er

!****

  elemental subroutine split_r (x, f, e)

    real(WP), intent(in)  :: x
    real(WP), intent(out) :: f
    integer, intent(out)  :: e

    ! Spit x into fraction and exponent parts

    f = FRACTION(x)
    e = EXPONENT(x)

    ! Finish

    return

  end subroutine split_r

end module gyre_ext_real
