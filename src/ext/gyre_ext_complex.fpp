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
    procedure             :: arith_${INFIX}_ec_ec
    procedure             :: arith_${INFIX}_ec_er
    procedure, pass(that) :: arith_${INFIX}_er_ec
    procedure             :: arith_${INFIX}_ec_r
    procedure, pass(that) :: arith_${INFIX}_r_ec
    procedure             :: arith_${INFIX}_ec_c
    procedure, pass(that) :: arith_${INFIX}_c_ec
    generic, public       :: operator($OP) => arith_${INFIX}_ec_ec, arith_${INFIX}_ec_er, arith_${INFIX}_er_ec, &
                                              arith_${INFIX}_ec_r, arith_${INFIX}_r_ec, &
                                              arith_${INFIX}_ec_c, arith_${INFIX}_c_ec
  $endsub

  $define $CMP_DECL $sub
    $local $INFIX $1
    $local $OP $2
    procedure             :: cmp_${INFIX}_ec_ec
    procedure             :: cmp_${INFIX}_ec_er
    procedure, pass(that) :: cmp_${INFIX}_er_ec
    procedure             :: cmp_${INFIX}_ec_r
    procedure, pass(that) :: cmp_${INFIX}_r_ec
    procedure             :: cmp_${INFIX}_ec_c
    procedure, pass(that) :: cmp_${INFIX}_c_ec
    generic, public       :: operator($OP) => cmp_${INFIX}_ec_ec, cmp_${INFIX}_ec_er, cmp_${INFIX}_er_ec, &
                                              cmp_${INFIX}_ec_r, cmp_${INFIX}_r_ec, &
                                              cmp_${INFIX}_ec_c, cmp_${INFIX}_c_ec
  $endsub

  type ext_complex_t
     private
     complex(WP) :: f ! Fractional part
     integer     :: e ! Exponent
   contains
     $ARITH_DECL(plus,+)
     procedure       :: arith_minus_ec
     generic, public :: operator(-) => arith_minus_ec
     $ARITH_DECL(minus,-)
     $ARITH_DECL(times,*)
     $ARITH_DECL(divide,/)
     $CMP_DECL(eq,==)
     $CMP_DECL(neq,/=)
  end type ext_complex_t

  ! Interface blocks

  $if($MPI)

  interface send
     module procedure send_ec_0
     module procedure send_ec_1
     module procedure send_ec_2
     module procedure send_ec_3
     module procedure send_ec_4
  end interface send

  interface recv
     module procedure recv_ec_0
     module procedure recv_ec_1
     module procedure recv_ec_2
     module procedure recv_ec_3
     module procedure recv_ec_4
  end interface recv

  interface recv_any
     module procedure recv_any_ec_0
     module procedure recv_any_ec_1
     module procedure recv_any_ec_2
     module procedure recv_any_ec_3
     module procedure recv_any_ec_4
  end interface recv_any

  interface bcast
     module procedure bcast_ec_0
     module procedure bcast_ec_1
     module procedure bcast_ec_2
     module procedure bcast_ec_3
     module procedure bcast_ec_4
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_ec_0
     module procedure bcast_alloc_ec_1
     module procedure bcast_alloc_ec_2
     module procedure bcast_alloc_ec_3
     module procedure bcast_alloc_ec_4
  end interface bcast_alloc

  interface gatherv
     module procedure gatherv_ec_0
     module procedure gatherv_ec_1
     module procedure gatherv_ec_2
     module procedure gatherv_ec_3
     module procedure gatherv_ec_4
  end interface gatherv

  interface allgatherv
     module procedure allgatherv_ec_0
     module procedure allgatherv_ec_1
     module procedure allgatherv_ec_2
     module procedure allgatherv_ec_3
     module procedure allgatherv_ec_4
  end interface allgatherv

  $endif

  interface ext_real
     module procedure ext_real_ec
  end interface ext_real

  interface ext_complex
     module procedure ext_complex_r
     module procedure ext_complex_er
     module procedure ext_complex_c
  end interface ext_complex

  interface real
     module procedure real_ec
  end interface real

  interface cmplx
     module procedure cmplx_ec
  end interface cmplx

  interface valid
     module procedure valid_ec
  end interface valid

  interface product
     module procedure product_ec
  end interface product

  interface exp
     module procedure exp_ec
  end interface exp

  interface abs
     module procedure abs_ec
  end interface abs

  interface fraction
     module procedure fraction_ec
  end interface fraction

  interface exponent
     module procedure exponent_ec
  end interface exponent

  interface scale
     module procedure scale_ec
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
  public :: ext_complex_t
  public :: ext_real
  public :: ext_complex
  public :: real
  public :: cmplx
  public :: valid
  public :: product
  public :: abs
  public :: exp
  public :: fraction
  public :: exponent
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

  $BCAST(ec,type(ext_complex_t),0)
  $BCAST(ec,type(ext_complex_t),1)
  $BCAST(ec,type(ext_complex_t),2)
  $BCAST(ec,type(ext_complex_t),3)
  $BCAST(ec,type(ext_complex_t),4)
  
!****

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

  $ALLGATHERV(ec,type(ext_complex_t),0)
  $ALLGATHERV(ec,type(ext_complex_t),1)
  $ALLGATHERV(ec,type(ext_complex_t),2)
  $ALLGATHERV(ec,type(ext_complex_t),3)
  $ALLGATHERV(ec,type(ext_complex_t),4)
  
  $endif

!****

  elemental function arith_plus_ec_ec (this, that) result (ez)

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
          f = this%f + that%f*RADIX_WP**(that%e - this%e)
          e = this%e
       else
          f = this%f*RADIX_WP**(this%e - that%e) + that%f
          e = that%e
       endif

       call split_c(f, ez%f, ez%e)
       ez%e = ez%e + e

    endif

    ! Finish

    return

  end function arith_plus_ec_ec

!****

  elemental function arith_minus_ec (this) result (ez)

    class(ext_complex_t), intent(in) :: this
    type(ext_complex_t)              :: ez

    ! Apply the unary minus operator

    ez%f = -this%f
    ez%e = this%e

    ! Finish

    return

  end function arith_minus_ec

!****

  elemental function arith_minus_ec_ec (this, that) result (ez)

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
          f = this%f - that%f*RADIX_WP**(that%e - this%e)
          e = this%e
       else
          f = this%f*RADIX_WP**(this%e - that%e) - that%f
          e = that%e
       endif

       call split_c(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_minus_ec_ec

!****

  elemental function arith_times_ec_ec (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the times operator

    if(this%f == 0._WP .OR. that%f == 0._WP) then

       ez = ext_complex(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_c(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_times_ec_ec

!****

  elemental function arith_divide_ec_ec (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    complex(WP) :: f
    integer     :: e

    ! Apply the divide operator

    if(this%f == 0._WP .AND. that%f /= 0._WP) then

       ez = ext_complex(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_c(f, ez%f, ez%e)
       ez = scale(ez, e)

    endif

    ! Finish

    return

  end function arith_divide_ec_ec

!****

  $define $ARITH $sub

  $local $INFIX $1
  $local $OP $2

  elemental function arith_${INFIX}_ec_er (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    class(ext_real_t), intent(in)    :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_er

  elemental function arith_${INFIX}_er_ec (this, that) result ($INFIX)

    class(ext_real_t), intent(in)    :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_er_ec

  elemental function arith_${INFIX}_ec_c (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    complex(WP), intent(in)          :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_c

  elemental function arith_${INFIX}_c_ec (this, that) result ($INFIX)

    complex(WP), intent(in)          :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_c_ec

  elemental function arith_${INFIX}_ec_r (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    real(WP), intent(in)             :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function arith_${INFIX}_ec_r

  elemental function arith_${INFIX}_r_ec (this, that) result ($INFIX)

    real(WP), intent(in)             :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: $INFIX

    ! Apply the arithmetic operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function arith_${INFIX}_r_ec

  $endsub

  $ARITH(plus,+)
  $ARITH(minus,-)
  $ARITH(times,*)
  $ARITH(divide,/)

!****

  elemental function cmp_eq_ec_ec (this, that) result (eq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: eq

    ! Apply the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function cmp_eq_ec_ec

!****

  elemental function cmp_neq_ec_ec (this, that) result (neq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: neq

    ! Apply the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function cmp_neq_ec_ec

!****

  $define $CMP $sub

  $local $INFIX $1
  $local $OP $2

  elemental function cmp_${INFIX}_ec_er (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    class(ext_real_t), intent(in)    :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_er

  elemental function cmp_${INFIX}_er_ec (this, that) result ($INFIX)

    class(ext_real_t), intent(in)    :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_er_ec

  elemental function cmp_${INFIX}_ec_r (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    real(WP), intent(in)             :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_r

  elemental function cmp_${INFIX}_r_ec (this, that) result ($INFIX)

    real(WP), intent(in)             :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_r_ec

  elemental function cmp_${INFIX}_ec_c (this, that) result ($INFIX)

    class(ext_complex_t), intent(in) :: this
    complex(WP), intent(in)          :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = this $OP ext_complex(that)

    ! Finish

    return

  end function cmp_${INFIX}_ec_c

  elemental function cmp_${INFIX}_c_ec (this, that) result ($INFIX)

    complex(WP), intent(in)          :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: $INFIX

    ! Apply the comparison operator

    $INFIX = ext_complex(this) $OP that

    ! Finish

    return

  end function cmp_${INFIX}_c_ec

  $endsub

  $CMP(eq,==)
  $CMP(neq,/=)

!****

  elemental function ext_real_ec (ez) result (ex)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)                :: ex

    ! Convert ext_complex to ext_real

    ex = ext_real(REAL(ez%f))
    ex = scale(ex, ez%e)

    ! Finish

    return

  end function ext_real_ec

!****

  elemental function ext_complex_r (x) result (ez)

    real(WP), intent(in) :: x
    type(ext_complex_t)  :: ez

    ! Convert real to ext_complex

    call split_c(CMPLX(x, KIND=WP), ez%f, ez%e)

    ! Finish

    return

  end function ext_complex_r

!****

  elemental function ext_complex_er (ex) result (ez)

    type(ext_real_t), intent(in) :: ex
    type(ext_complex_t)          :: ez

    ! Convert ext_real to ext_complex

    call split_c(CMPLX(FRACTION(ex), KIND=WP), ez%f, ez%e)
    ez = scale(ez, EXPONENT(ex))

    ! Finish

    return

  end function ext_complex_er

!****

  elemental function ext_complex_c (z) result (ez)

    complex(WP), intent(in) :: z
    type(ext_complex_t)     :: ez

    ! Convert complex to ext_complex

    call split_c(z, ez%f, ez%e)

    ! Finish

    return

  end function ext_complex_c

!****
    
  elemental function real_ec (ez) result (x)

    type(ext_complex_t), intent(in) :: ez
    real(WP)                        :: x

    ! Convert ext_complex to real

    x = REAL(cmplx(ez))

    ! Finish

    return

  end function real_ec

!****
    
  elemental function cmplx_ec (ez) result (z)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: z

    ! Convert ext_complex to complex

    if(ez%f /= 0._WP) then
       z = ez%f*RADIX_WP**ez%e
    else
       z = 0._WP
    endif

    ! Finish

    return

  end function cmplx_ec

!****

  elemental function valid_ec (ez) result (valid)

    type(ext_complex_t), intent(in) :: ez
    logical                         :: valid

    ! Determine if ext_complex is valid

    valid = ABS(REAL(ez%f)) >= RADIX_WP**(-1) .AND. ABS(REAL(ez%f)) < 1._WP

    ! Finish

    return

  end function valid_ec

!****

  function product_ec (ez) result (prod_ez)

    type(ext_complex_t), intent(in) :: ez(:)
    type(ext_complex_t)             :: prod_ez

    integer     :: i
    complex(WP) :: f
    integer     :: e
    
    ! Calculate the product of the elements of ez

    prod_ez%f = 1._WP
    prod_ez%e = SUM(ez%e)

    do i = 1,SIZE(ez)
       
       call split_c(prod_ez%f*ez(i)%f, f, e)

       prod_ez%f = f
       prod_ez = scale(prod_ez, e)

    end do

    ! Finish

    return

  end function product_ec

!****

  elemental function abs_ec (ez) result (abs_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)               :: abs_ez

    ! Calculate the absolute value of ez

    abs_ez = ext_real(ABS(ez%f))
    abs_ez = scale(abs_ez, ez%e)

    ! Finish

    return

  end function abs_ec

!****

  elemental function exp_ec (ez) result (exp_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_complex_t)             :: exp_ez

    type(ext_real_t) :: exp_ex

    ! Calculate the exponential of ez

    exp_ex = exp(ext_real(ez))

    exp_ez%f = FRACTION(exp_ex)*EXP((0._WP,1._WP)*AIMAG(cmplx(ez)))
    exp_ez%e = EXPONENT(exp_ex)

    ! Finish

    return

  end function exp_ec

!****

  elemental function fraction_ec (ez) result (fraction_ez)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: fraction_ez

    ! Return the fraction part of ez

    fraction_ez = ez%f

    ! Finish

    return

  end function fraction_ec

!****

  elemental function exponent_ec (ez) result (exponent_ez)

    type(ext_complex_t), intent(in) :: ez
    integer                         :: exponent_ez

    ! Return the exponent part of ez

    exponent_ez = ez%e

    ! Finish

    return

  end function exponent_ec

!****

  elemental function scale_ec (ez, de) result (scale_ez)

    class(ext_complex_t), intent(in) :: ez
    integer, intent(in)              :: de
    type(ext_complex_t)              :: scale_ez

    ! Scale ez by RADIX_WP**de

    scale_ez%f = ez%f
    scale_ez%e = ez%e + de

    ! Finish

    return

  end function scale_ec

!****

  elemental subroutine split_c (z, f, e)

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

  end subroutine split_c

end module gyre_ext_complex
