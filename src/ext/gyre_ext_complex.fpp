! Module   : gyre_ezt_complex
! Purpose  : extented-range arithmetic (complex)
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

module gyre_ext_complex

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_ext_real

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  real(WP), parameter :: RADIX_WP = REAL(RADIX(1._WP), WP)
  
  ! Derived-type definitions

  $define $ARITH_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: arith_${INFIX}_ec_ec_
    procedure             :: arith_${INFIX}_ec_er_
    procedure, pass(that) :: arith_${INFIX}_er_ec_
    procedure             :: arith_${INFIX}_ec_r_
    procedure, pass(that) :: arith_${INFIX}_r_ec_
    procedure             :: arith_${INFIX}_ec_c_
    procedure, pass(that) :: arith_${INFIX}_c_ec_
    generic, public       :: operator($OP) => arith_${INFIX}_ec_ec_, arith_${INFIX}_ec_er_, arith_${INFIX}_er_ec_, &
                                              arith_${INFIX}_ec_r_, arith_${INFIX}_r_ec_, &
                                              arith_${INFIX}_ec_c_, arith_${INFIX}_c_ec_
  $endsub

  $define $CMP_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: cmp_${INFIX}_ec_ec_
    procedure             :: cmp_${INFIX}_ec_er_
    procedure, pass(that) :: cmp_${INFIX}_er_ec_
    procedure             :: cmp_${INFIX}_ec_r_
    procedure, pass(that) :: cmp_${INFIX}_r_ec_
    procedure             :: cmp_${INFIX}_ec_c_
    procedure, pass(that) :: cmp_${INFIX}_c_ec_
    generic, public       :: operator($OP) => cmp_${INFIX}_ec_ec_, cmp_${INFIX}_ec_er_, cmp_${INFIX}_er_ec_, &
                                              cmp_${INFIX}_ec_r_, cmp_${INFIX}_r_ec_, &
                                              cmp_${INFIX}_ec_c_, cmp_${INFIX}_c_ec_
  $endsub

  type ext_complex_t
     private
     complex(WP) :: f ! Fractional part
     integer     :: e ! Exponent
   contains
     $ARITH_DECL(plus,+)
     procedure       :: arith_minus_ec_
     generic, public :: operator(-) => arith_minus_ec_
     $ARITH_DECL(minus,-)
     $ARITH_DECL(times,*)
     $ARITH_DECL(divide,/)
     $CMP_DECL(eq,==)
     $CMP_DECL(neq,/=)
  end type ext_complex_t

  ! Interface blocks

  interface ext_complex_t
     module procedure ext_complex_t_r_
     module procedure ext_complex_t_c_
     module procedure ext_complex_t_er_
  end interface ext_complex_t

  interface ext_real_t
     module procedure ext_real_t_ec_
  end interface ext_real_t

  interface real
     module procedure real_
  end interface real

  interface cmplx
     module procedure cmplx_
  end interface cmplx

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

  public :: ext_complex_t
  public :: ext_real_t
  public :: real
  public :: cmplx
  public :: valid
  public :: product
  public :: abs
  public :: exp
  public :: sqrt
  public :: fraction
  public :: exponent
  public :: scale
  $if($MPI)
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

  elemental function ext_complex_t_r_ (x) result (ez)

    real(WP), intent(in) :: x
    type(ext_complex_t)  :: ez

    ! Construct the ext_complex_t from the real x

    call split_(CMPLX(x, KIND=WP), ez%f, ez%e)

    ! Finish

    return

  end function ext_complex_t_r_

!****

  elemental function ext_complex_t_c_ (z) result (ez)

    complex(WP), intent(in) :: z
    type(ext_complex_t)     :: ez

    ! Construct the ext_complex_t from the complex z

    call split_(z, ez%f, ez%e)

    ! Finish

    return

  end function ext_complex_t_c_

!****

  elemental function ext_complex_t_er_ (ex) result (ez)

    type(ext_real_t), intent(in) :: ex
    type(ext_complex_t)          :: ez

    ! Construct the ext_complex_t from the ext_real_t ex

    call split_(CMPLX(FRACTION(ex), KIND=WP), ez%f, ez%e)
    ez = scale(ez, EXPONENT(ex))

    ! Finish

    return

  end function ext_complex_t_er_

!****

  elemental function arith_plus_ec_ec_ (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the plus operator

    if(this%f == 0._WP) then

       ez%f = that%f
       ez%e = that%e

    elseif(that%f == 0._WP) then

       ez%f = this%f
       ez%e = this%e

    else

       if(this%e > that%e) then
          f = this%f + cmplx(ext_complex_t(that%f, that%e - this%e))
          e = this%e
       else
          f = cmplx(ext_complex_t(this%f, this%e - that%e)) + that%f
          e = that%e
       endif

       call split_(f, ez%f, ez%e)
       ez%e = ez%e + e

    endif

    ! Finish

    return

  end function arith_plus_ec_ec_

!****

  elemental function arith_minus_ec_ (this) result (ez)

    class(ext_complex_t), intent(in) :: this
    type(ext_complex_t)              :: ez

    ! Apply the unary minus operator

    ez%f = -this%f
    ez%e = this%e

    ! Finish

    return

  end function arith_minus_ec_

!****

  elemental function arith_minus_ec_ec_ (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the minus operator

    if(this%f == 0._WP) then
       
       ez%f = -that%f
       ez%e = that%e

    elseif(that%f == 0._WP) then

       ez%f = this%f
       ez%e = this%e

    else

       if(this%e > that%e) then
          f = this%f - cmplx(ext_complex_t(that%f, that%e - this%e))
          e = this%e
       else
          f = cmplx(ext_complex_t(this%f, this%e - that%e)) - that%f
          e = that%e
       endif

       call split_(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_minus_ec_ec_

!****

  elemental function arith_times_ec_ec_ (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the times operator

    if(this%f == 0._WP .OR. that%f == 0._WP) then

       ez = ext_complex_t(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_times_ec_ec_

!****

  elemental function arith_divide_ec_ec_ (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the divide operator

    if(this%f == 0._WP .AND. that%f /= 0._WP) then

       ez = ext_complex_t(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_divide_ec_ec_

!****

  $define $ARITH $sub

  $local $INFIX $1
  $local $OP $2

  elemental function arith_${INFIX}_ec_er_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    class(ext_real_t), intent(in)    :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_er_

  elemental function arith_${INFIX}_er_ec_ (this, that) result ($INFIX)

    class(ext_real_t), intent(in)    :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_er_ec_

  elemental function arith_${INFIX}_ec_c_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    complex(WP), intent(in)          :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_c_

  elemental function arith_${INFIX}_c_ec_ (this, that) result ($INFIX)

    complex(WP), intent(in)          :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_c_ec_

  elemental function arith_${INFIX}_ec_r_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    real(WP), intent(in)             :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_r_

  elemental function arith_${INFIX}_r_ec_ (this, that) result ($INFIX)

    real(WP), intent(in)             :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_r_ec_

  $endsub

  $ARITH(plus,+)
  $ARITH(minus,-)
  $ARITH(times,*)
  $ARITH(divide,/)

!****

  elemental function cmp_eq_ec_ec_ (this, that) result (eq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: eq

    ! Apply the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function cmp_eq_ec_ec_

!****

  elemental function cmp_neq_ec_ec_ (this, that) result (neq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: neq

    ! Apply the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function cmp_neq_ec_ec_

!****

  $define $CMP $sub

  $local $INFIX $1
  $local $OP $2

  elemental function cmp_${INFIX}_ec_er_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    class(ext_real_t), intent(in)    :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_er_

  elemental function cmp_${INFIX}_er_ec_ (this, that) result ($INFIX)

    class(ext_real_t), intent(in)    :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_er_ec_

  elemental function cmp_${INFIX}_ec_r_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    real(WP), intent(in)             :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_r_

  elemental function cmp_${INFIX}_r_ec_ (this, that) result ($INFIX)

    real(WP), intent(in)             :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_r_ec_

  elemental function cmp_${INFIX}_ec_c_ (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    complex(WP), intent(in)          :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex_t(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_c_

  elemental function cmp_${INFIX}_c_ec_ (this, that) result ($INFIX)

    complex(WP), intent(in)          :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex_t(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_c_ec_

  $endsub

  $CMP(eq,==)
  $CMP(neq,/=)

!****

  elemental function ext_real_t_ec_ (ez) result (ex)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)                :: ex

    ! Construct the ext_real_t from the ext_complex_t ez

    ex = ext_real_t(REAL(ez%f))
    ex = scale(ex, ez%e)

    ! Finish

    return

  end function ext_real_t_ec_

!****
    
  elemental function real_ (ez) result (x)

    type(ext_complex_t), intent(in) :: ez
    real(WP)                        :: x

    ! Convert ext_complex_t to real

    x = REAL(cmplx(ez))

    ! Finish

    return

  end function real_

!****
    
  elemental function cmplx_ (ez) result (z)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: z

    integer :: e_min

    ! Convert ext_complex_t to complex

    if(ez%f /= 0._WP) then

       e_min = MINEXPONENT(0._WP)

       if(ez%e >= e_min) then
          z = ez%f*RADIX_WP**ez%e
       else
          z = (ez%f*RADIX_WP**MAX(ez%e-e_min, -DIGITS(0._WP)-1))*RADIX_WP**e_min
       endif

    else

       z = 0._WP

    endif

    ! Finish

    return

  end function cmplx_

!****

  elemental function valid_ (ez) result (valid)

    type(ext_complex_t), intent(in) :: ez
    logical                         :: valid

    ! Determine if ez is valid

    valid = ABS(REAL(ez%f)) >= RADIX_WP**(-1) .AND. ABS(REAL(ez%f)) < 1._WP

    ! Finish

    return

  end function valid_

!****

  function product_ (ez) result (prod_ez)

    type(ext_complex_t), intent(in) :: ez(:)
    type(ext_complex_t)             :: prod_ez

    integer     :: i
    complex(WP) :: f
    integer     :: e
    
    ! Calculate the product of the elements of ez

    prod_ez%f = 1._WP
    prod_ez%e = SUM(ez%e)

    do i = 1,SIZE(ez)
       
       call split_(prod_ez%f*ez(i)%f, f, e)

       prod_ez%f = f
       prod_ez = scale(prod_ez, e)

    end do

    ! Finish

    return

  end function product_

!****

  elemental function abs_ (ez) result (abs_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)               :: abs_ez

    ! Calculate the absolute value of ez

    abs_ez = ext_real_t(ABS(ez%f))
    abs_ez = scale(abs_ez, ez%e)

    ! Finish

    return

  end function abs_

!****

  elemental function exp_ (ez) result (exp_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_complex_t)             :: exp_ez

    type(ext_real_t) :: exp_ex

    ! Calculate the exponential of ez

    exp_ex = exp(ext_real_t(ez))

    exp_ez%f = FRACTION(exp_ex)*EXP((0._WP,1._WP)*AIMAG(cmplx(ez)))
    exp_ez%e = EXPONENT(exp_ex)

    ! Finish

    return

  end function exp_

!****

  elemental function sqrt_ (ez) result (sqrt_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_complex_t)             :: sqrt_ez

    ! Calculate the square root of ez

    sqrt_ez = ext_complex_t(SQRT(ez%f))

    if (MOD(ez%e, 2) == 0) then
       sqrt_ez = scale(sqrt_ez, ez%e/2)
    else
       sqrt_ez = scale(sqrt_ez, (ez%e-1)/2)*SQRT(2._WP)
    endif

    ! Finish

    return

  end function sqrt_

!****

  elemental function fraction_ (ez) result (fraction_ez)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: fraction_ez

    ! Return the fraction part of ez

    fraction_ez = ez%f

    ! Finish

    return

  end function fraction_

!****

  elemental function exponent_ (ez) result (exponent_ez)

    type(ext_complex_t), intent(in) :: ez
    integer                         :: exponent_ez

    ! Return the exponent part of ez

    exponent_ez = ez%e

    ! Finish

    return

  end function exponent_

!****

  elemental function scale_ (ez, de) result (scale_ez)

    class(ext_complex_t), intent(in) :: ez
    integer, intent(in)              :: de
    type(ext_complex_t)              :: scale_ez

    ! Scale ez by RADIX_WP**de

    scale_ez%f = ez%f
    scale_ez%e = ez%e + de

    ! Finish

    return

  end function scale_

!****

  elemental subroutine split_ (z, f, e)

    complex(WP), intent(in)  :: z
    complex(WP), intent(out) :: f
    integer, intent(out)     :: e

    real(WP)    :: z_r
    real(WP)    :: z_i
    real(WP)    :: f_r
    real(WP)    :: f_i
    integer     :: e_r
    integer     :: e_i
    integer     :: e_ref

    ! Spit z into fraction and exponent parts

    z_r = REAL(z)
    z_i = AIMAG(z)

    f_r = FRACTION(z_r)
    f_i = FRACTION(z_i)

    if(f_r == 0._WP .AND. f_i == 0._WP) then

       f = 0._WP
       e = 0

    else

       e_r = EXPONENT(z_r)
       e_i = EXPONENT(z_i)

       e_ref = MAXVAL([e_r,e_i], MASK=[f_r,f_i] /= 0._WP)

       if(e_r - e_ref < MINEXPONENT(0._WP)) then
          f_r = 0
       else
          f_r = f_r*RADIX_WP**(e_r - e_ref)
       endif
       
       if(e_i - e_ref < MINEXPONENT(0._WP)) then
          f_i = 0
       else
          f_i = f_i*RADIX_WP**(e_i - e_ref)
       endif

       f = CMPLX(f_r, f_i, WP)
       e = e_ref

    endif

    ! Finish

    return

  end subroutine split_

!****

  $if($MPI)

  $define $SEND $sub

  $local $BUFFER_RANK $1

  subroutine send_${BUFFER_RANK}_ (buffer, dest_rank, tag, sync)

    type(ext_complex_t), intent(in) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)             :: dest_rank
    integer, optional, intent(in)   :: tag
    logical, optional, intent(in)   :: sync

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

    type(ext_complex_t), intent(out) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)              :: src_rank
    integer, optional, intent(in)    :: tag

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

    type(ext_complex_t), intent(out) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(out)             :: src_rank
    integer, optional, intent(in)    :: tag

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

    type(ext_complex_t), intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)                :: root_rank

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

  $BCAST_ALLOC(type(ext_complex_t),0)
  $BCAST_ALLOC(type(ext_complex_t),1)
  $BCAST_ALLOC(type(ext_complex_t),2)
  $BCAST_ALLOC(type(ext_complex_t),3)
  $BCAST_ALLOC(type(ext_complex_t),4)

!****

  $define $GATHERV $sub

  $local $BUFFER_RANK $1

  subroutine gatherv_${BUFFER_RANK}_ (send_buffer, sendcount, recv_buffer, recvcounts, displs, root_rank)

    type(ext_complex_t), intent(in)    :: send_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)                :: sendcount
    type(ext_complex_t), intent(inout) :: recv_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)                :: recvcounts(:)
    integer, intent(in)                :: displs(:)
    integer, intent(in)                :: root_rank

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

    type(ext_complex_t), intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)                :: recvcounts(:)
    integer, intent(in)                :: displs(:)

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

end module gyre_ext_complex
