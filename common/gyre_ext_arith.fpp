! Module   : gyre_ext_arith
! Purpose  : extented-range arithmetic
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

module gyre_ext_arith

  ! Uses

  use core_kinds
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  real(WP), parameter :: RADIX_WP = REAL(RADIX(1._WP), WP)
  
  ! Derived-type definitions

  type ext_real_t
     private
     real(WP) :: f ! Fractional part
     integer  :: e ! Exponent
   contains
     private
     procedure       :: unary_minus_er
     procedure       :: binary_times_er
     procedure       :: binary_eq_er
     procedure       :: binary_neq_er
     procedure       :: binary_lt_er
     procedure       :: binary_gt_er
     procedure       :: binary_le_er
     procedure       :: binary_ge_er
     generic, public :: operator(-) => unary_minus_er
     generic, public :: operator(*) => binary_times_er
     generic, public :: operator(==) => binary_eq_er
     generic, public :: operator(/=) => binary_neq_er
     generic, public :: operator(<) => binary_lt_er
     generic, public :: operator(>) => binary_gt_er
     generic, public :: operator(<=) => binary_le_er
     generic, public :: operator(>=) => binary_ge_er
  end type ext_real_t

  type ext_complex_t
     private
     complex(WP) :: f ! Fractional part
     integer     :: e ! Exponent
   contains
     private
     procedure       :: unary_minus_ec
     procedure       :: binary_times_ec
     procedure       :: binary_eq_ec
     procedure       :: binary_neq_ec
     generic, public :: operator(-) => unary_minus_ec
     generic, public :: operator(*) => binary_times_ec
     generic, public :: operator(==) => binary_eq_ec
     generic, public :: operator(/=) => binary_neq_ec
  end type ext_complex_t

  ! Interface blocks

  $if($MPI)

  interface send
     module procedure send_er_0
     module procedure send_er_1
     module procedure send_er_2
     module procedure send_er_3
     module procedure send_er_4
     module procedure send_ec_0
     module procedure send_ec_1
     module procedure send_ec_2
     module procedure send_ec_3
     module procedure send_ec_4
  end interface send

  interface recv
     module procedure recv_er_0
     module procedure recv_er_1
     module procedure recv_er_2
     module procedure recv_er_3
     module procedure recv_er_4
     module procedure recv_ec_0
     module procedure recv_ec_1
     module procedure recv_ec_2
     module procedure recv_ec_3
     module procedure recv_ec_4
  end interface recv

  interface recv_any
     module procedure recv_any_er_0
     module procedure recv_any_er_1
     module procedure recv_any_er_2
     module procedure recv_any_er_3
     module procedure recv_any_er_4
     module procedure recv_any_ec_0
     module procedure recv_any_ec_1
     module procedure recv_any_ec_2
     module procedure recv_any_ec_3
     module procedure recv_any_ec_4
  end interface recv_any

  interface bcast
     module procedure bcast_er_0
     module procedure bcast_er_1
     module procedure bcast_er_2
     module procedure bcast_er_3
     module procedure bcast_er_4
     module procedure bcast_ec_0
     module procedure bcast_ec_1
     module procedure bcast_ec_2
     module procedure bcast_ec_3
     module procedure bcast_ec_4
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_er_0
     module procedure bcast_alloc_er_1
     module procedure bcast_alloc_er_2
     module procedure bcast_alloc_er_3
     module procedure bcast_alloc_er_4
     module procedure bcast_alloc_ec_0
     module procedure bcast_alloc_ec_1
     module procedure bcast_alloc_ec_2
     module procedure bcast_alloc_ec_3
     module procedure bcast_alloc_ec_4
  end interface bcast_alloc

  interface gatherv
     module procedure gatherv_er_0
     module procedure gatherv_er_1
     module procedure gatherv_er_2
     module procedure gatherv_er_3
     module procedure gatherv_er_4
     module procedure gatherv_ec_0
     module procedure gatherv_ec_1
     module procedure gatherv_ec_2
     module procedure gatherv_ec_3
     module procedure gatherv_ec_4
  end interface gatherv

  interface allgatherv
     module procedure allgatherv_er_0
     module procedure allgatherv_er_1
     module procedure allgatherv_er_2
     module procedure allgatherv_er_3
     module procedure allgatherv_er_4
     module procedure allgatherv_ec_0
     module procedure allgatherv_ec_1
     module procedure allgatherv_ec_2
     module procedure allgatherv_ec_3
     module procedure allgatherv_ec_4
  end interface allgatherv

  $endif

  interface ext_real
     module procedure ext_real_r
     module procedure ext_real_c
     module procedure ext_real_ec
  end interface ext_real

  interface ext_complex
     module procedure ext_complex_r
     module procedure ext_complex_er
     module procedure ext_complex_c
  end interface ext_complex

  interface real
     module procedure real_r
     module procedure real_c
  end interface real

  interface cmplx
     module procedure cmplx_r
     module procedure cmplx_c
  end interface cmplx

  interface valid
     module procedure valid_r
     module procedure valid_c
  end interface valid

  interface product
     module procedure product_r
     module procedure product_c
  end interface product

  interface exp
     module procedure exp_r
     module procedure exp_c
  end interface exp

  interface abs
     module procedure abs_r
     module procedure abs_c
  end interface abs

  interface fraction
     module procedure fraction_r
     module procedure fraction_c
  end interface fraction

  interface exponent
     module procedure exponent_r
     module procedure exponent_c
  end interface exponent

  interface scale
     module procedure scale_r
     module procedure scale_c
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

  elemental function ext_real_ec (ez) result (ex)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)                :: ex

    ! Convert ext_complex to ext_real

    call split_r(REAL(ez%f), ex%f, ex%e)
    ex%e = ex%e + ez%e

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

    call split_c(CMPLX(ex%f, KIND=WP), ez%f, ez%e)
    ez%e = ez%e + ex%e

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

  elemental function real_r (ex) result (x)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: x

    ! Convert ext_real to real

    x = ex%f*RADIX_WP**ex%e

    ! Finish

    return

  end function real_r

!****
    
  elemental function real_c (ez) result (x)

    type(ext_complex_t), intent(in) :: ez
    real(WP)                        :: x

    ! Convert ext_complex to real

    x = REAL(cmplx(ez))

    ! Finish

    return

  end function real_c

!****
    
  elemental function cmplx_r (ex) result (z)

    type(ext_real_t), intent(in) :: ex
    complex(WP)                  :: z

    ! Convert ext_real to complex

    z = CMPLX(real(ex), KIND=WP)

    ! Finish

    return

  end function cmplx_r

!****
    
  elemental function cmplx_c (ez) result (z)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: z

    ! Convert ext_complex to complex

    z = CMPLX(ez%f, KIND=WP)*RADIX_WP**ez%e

    ! Finish

    return

  end function cmplx_c

!****

  elemental function valid_r (ex) result (valid)

    type(ext_real_t), intent(in) :: ex
    logical                      :: valid

    ! Determine if ext_real is valid

    valid = ABS(ex%f) >= RADIX_WP**(-1) .AND. ABS(ex%f) < 1._WP

    ! Finish

    return

  end function valid_r

!****

  elemental function valid_c (ez) result (valid)

    type(ext_complex_t), intent(in) :: ez
    logical                         :: valid

    ! Determine if ext_complex is valid

    valid = ABS(REAL(ez%f)) >= RADIX_WP**(-1) .AND. ABS(REAL(ez%f)) < 1._WP

    ! Finish

    return

  end function valid_c

!****

  function product_r (ex) result (prod_ex)

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
       prod_ex%e = prod_ex%e + EXPONENT(g)

    end do

    ! Finish

    return

  end function product_r

!****

  function product_c (ez) result (prod_ez)

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
       prod_ez%e = prod_ez%e + e

    end do

    ! Finish

    return

  end function product_c

!****

  elemental function abs_r (ex) result (abs_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: abs_ex

    ! Calculate the absolute value of ex

    abs_ex%f = ABS(ex%f)
    abs_ex%e = ex%e

    ! Finish

    return

  end function abs_r

!****

  elemental function abs_c (ez) result (abs_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_real_t)                :: abs_ez

    ! Calculate the absolute value of ez

    abs_ez%f = ABS(ez%f)
    abs_ez%e = ez%e

    ! Finish

    return

  end function abs_c

!****

  elemental function exp_r (ex) result (exp_ex)

    type(ext_real_t), intent(in) :: ex
    type(ext_real_t)             :: exp_ex

    real(WP) :: g
    integer  :: e

    ! Calculate the exponential of ex

    g = real(ex)/LOG(RADIX_WP)
    e = FLOOR(g)

    call split_r(RADIX_WP**(g-e), exp_ex%f, exp_ex%e)

    exp_ex%e = exp_ex%e + e

    ! Finish

    return

  end function exp_r

!****

  elemental function exp_c (ez) result (exp_ez)

    type(ext_complex_t), intent(in) :: ez
    type(ext_complex_t)             :: exp_ez

    type(ext_real_t) :: exp_ex

    ! Calculate the exponential of ez

    exp_ex = exp_r(ext_real(ez))

    exp_ez%f = exp_ex%f*EXP((0._WP,1._WP)*AIMAG(cmplx(ez)))
    exp_ez%e = exp_ex%e

    ! Finish

    return

  end function exp_c

!****

  elemental function fraction_r (ex) result (fraction_ex)

    type(ext_real_t), intent(in) :: ex
    real(WP)                     :: fraction_ex

    ! Return the fraction part of ex

    fraction_ex = ex%f

    ! Finish

    return

  end function fraction_r

!****

  elemental function fraction_c (ez) result (fraction_ez)

    type(ext_complex_t), intent(in) :: ez
    complex(WP)                     :: fraction_ez

    ! Return the fraction part of ez

    fraction_ez = ez%f

    ! Finish

    return

  end function fraction_c

!****

  elemental function exponent_r (ex) result (exponent_ex)

    type(ext_real_t), intent(in) :: ex
    integer                      :: exponent_ex

    ! Return the exponent part of ex

    exponent_ex = ex%e

    ! Finish

    return

  end function exponent_r

!****

  elemental function exponent_c (ez) result (exponent_ez)

    type(ext_complex_t), intent(in) :: ez
    integer                         :: exponent_ez

    ! Return the exponent part of ez

    exponent_ez = ez%e

    ! Finish

    return

  end function exponent_c

!****

  elemental function scale_r (ex, de) result (scale_ex)

    class(ext_real_t), intent(in) :: ex
    integer, intent(in)           :: de
    type(ext_real_t)              :: scale_ex

    ! Scale ex by 2**de

    scale_ex%f = ex%f
    scale_ex%e = ex%e + de

    ! Finish

    return

  end function scale_r

!****

  elemental function scale_c (ez, de) result (scale_ez)

    class(ext_complex_t), intent(in) :: ez
    integer, intent(in)              :: de
    type(ext_complex_t)              :: scale_ez

    ! Scale ez by 2**de

    scale_ez%f = ez%f
    scale_ez%e = ez%e + de

    ! Finish

    return

  end function scale_c

!****

  elemental function unary_minus_er (this) result (ex)

    class(ext_real_t), intent(in) :: this
    type(ext_real_t)              :: ex

    ! Apply the unary minus operator

    ex%f = -this%f
    ex%e = this%e

    ! Finish

    return

  end function unary_minus_er

!****

  elemental function binary_times_er (this, that) result (ex)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    type(ext_real_t)              :: ex

    ! Apply the binary times operator

    call split_r(this%f*that%f, ex%f, ex%e)
    ex%e = ex%e + this%e + that%e

    ! Finish

    return

  end function binary_times_er

!****

  elemental function binary_eq_er (this, that) result (eq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: eq

    ! Apply the binary equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function binary_eq_er

!****

  elemental function binary_neq_er (this, that) result (neq)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: neq

    ! Apply the binary inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function binary_neq_er

!****

  elemental function binary_lt_er (this, that) result (lt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: lt

    real(WP) :: this_s
    real(WP) :: that_s

    ! Apply the binary less-than operator

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

  end function binary_lt_er

!****

  elemental function binary_gt_er (this, that) result (gt)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: gt

    ! Apply the binary greater-than operator

    gt = that < this

    ! Finish

    return

  end function binary_gt_er

!****

  elemental function binary_le_er (this, that) result (le)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: le

    ! Apply the binary less-than-or-equal operator

    le = .NOT. this > that

    ! Finish

    return

  end function binary_le_er
  
!****

  elemental function binary_ge_er (this, that) result (ge)

    class(ext_real_t), intent(in) :: this
    class(ext_real_t), intent(in) :: that
    logical                       :: ge

    ! Apply the binary greater-than-or-equal operator

    ge = .NOT. this < that

    ! Finish

    return

  end function binary_ge_er

!****

  elemental function unary_minus_ec (this) result (ez)

    class(ext_complex_t), intent(in) :: this
    type(ext_complex_t)              :: ez

    ! Apply the unary minus operator

    ez%f = -this%f
    ez%e = this%e

    ! Finish

    return

  end function unary_minus_ec

!****

  elemental function binary_times_ec (this, that) result (ez)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    type(ext_complex_t)              :: ez

    ! Apply the binary times operator

    call split_c(this%f*that%f, ez%f, ez%e)
    ez%e = ez%e + this%e + that%e

    ! Finish

    return

  end function binary_times_ec

!****

  elemental function binary_eq_ec (this, that) result (eq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: eq
 
    ! Apply the binary equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function binary_eq_ec

!****

  elemental function binary_neq_ec (this, that) result (neq)

    class(ext_complex_t), intent(in) :: this
    class(ext_complex_t), intent(in) :: that
    logical                          :: neq

    ! Apply the binary inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function binary_neq_ec

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

end module gyre_ext_arith
