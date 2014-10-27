! Module   : gyre_sysmtx
! Purpose  : blocked system matrix
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

module gyre_sysmtx

  ! Uses

  use core_kinds
  use core_parallel
  use core_linalg

  use gyre_ext_arith
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type sysmtx_t
     private
     complex(WP), allocatable         :: B_i(:,:)   ! Inner boundary conditions
     complex(WP), allocatable         :: B_o(:,:)   ! Outer boundary conditions
     complex(WP), allocatable         :: E_l(:,:,:) ! Left equation blocks
     complex(WP), allocatable         :: E_r(:,:,:) ! Right equation blocks
     type(ext_complex_t)              :: S_i        ! Inner boundary scale
     type(ext_complex_t)              :: S_o        ! Outer boundary scale
     type(ext_complex_t), allocatable :: S(:)       ! Block scales
     integer                          :: n          ! Number of equation blocks
     integer                          :: n_e        ! Number of equations per block
     integer                          :: n_i        ! Number of inner boundary conditions
     integer                          :: n_o        ! Number of outer boundary conditions
   contains
     private
     procedure, public :: set_inner_bound => set_inner_bound_
     procedure, public :: set_outer_bound => set_outer_bound_
     procedure, public :: set_block => set_block_
     procedure, public :: scale_rows => scale_rows_
     procedure, public :: determinant => determinant_
     procedure, public :: null_vector => null_vector_
     procedure, public :: sol_vector => sol_vector_
  end type sysmtx_t

  ! Interfaces

  interface sysmtx_t
     module procedure sysmtx_t_
  end interface sysmtx_t

  ! Access specifiers

  private

  public :: sysmtx_t

  ! Procedures

contains

  function sysmtx_t_ (n, n_e, n_i, n_o) result (sm)

    integer, intent(in) :: n
    integer, intent(in) :: n_e
    integer, intent(in) :: n_i
    integer, intent(in) :: n_o
    type(sysmtx_t)      :: sm

    ! Construct the sysmtx_t

    allocate(sm%E_l(n_e,n_e,n))
    allocate(sm%E_r(n_e,n_e,n))

    allocate(sm%B_i(n_i,n_e))
    allocate(sm%B_o(n_o,n_e))

    allocate(sm%S(n))

    sm%n = n
    sm%n_e = n_e
    sm%n_i = n_i
    sm%n_o = n_o

    ! Finish

    return

  end function sysmtx_t_

!****

  subroutine set_inner_bound_ (this, B_i, S_i)

    class(sysmtx_t), intent(inout)  :: this
    complex(WP), intent(in)         :: B_i(:,:)
    type(ext_complex_t), intent(in) :: S_i
    
    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    ! Set the inner boundary conditions

    this%B_i = B_i
    this%S_i = S_i

    ! Finish

    return

  end subroutine set_inner_bound_

!****

  subroutine set_outer_bound_ (this, B_o, S_o)

    class(sysmtx_t), intent(inout)  :: this
    complex(WP), intent(in)         :: B_o(:,:)
    type(ext_complex_t), intent(in) :: S_o

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Set the outer boundary conditions

    this%B_o = B_o
    this%S_o = S_o

    ! Finish

    return

  end subroutine set_outer_bound_

!****

  subroutine set_block_ (this, k, E_l, E_r, S)

    class(sysmtx_t), intent(inout)  :: this
    integer, intent(in)             :: k
    complex(WP), intent(in)         :: E_l(:,:)
    complex(WP), intent(in)         :: E_r(:,:)
    type(ext_complex_t), intent(in) :: S

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    $ASSERT(k >= 1,Invalid block index)
    $ASSERT(k <= this%n,Invalid block index)

    ! Set the block

    this%E_l(:,:,k) = E_l
    this%E_r(:,:,k) = E_r

    this%S(k) = S

    ! Finish

    return

  end subroutine set_block_

!****

  subroutine scale_rows_ (this)

    class(sysmtx_t), intent(inout) :: this

    real(WP) :: scale
    integer  :: i
    integer  :: k

    ! Scale the rows of the sysmtx to have maximum absolute value of unity

    do i = 1,this%n_i
       scale = MAXVAL(ABS(this%B_i(i,:)))
       this%B_i(i,:) = this%B_i(i,:)/scale
       this%S_i = this%S_i*scale
    end do

    do k = 1, this%n
       do i = 1, this%n_e
          scale = MAX(MAXVAL(ABS(this%E_l(i,:,k))), MAXVAL(ABS(this%E_r(i,:,k))))
          this%E_l(i,:,k) = this%E_l(i,:,k)/scale
          this%E_r(i,:,k) = this%E_r(i,:,k)/scale
          this%S(k) = this%S(k)*scale
       end do
    end do

    do i = 1,this%n_o
       scale = MAXVAL(ABS(this%B_o(i,:)))
       this%B_o(i,:) = this%B_o(i,:)/scale
       this%S_o = this%S_o*scale
    end do

    ! Finish

    return

  end subroutine scale_rows_

!****

  subroutine determinant_ (this, det, use_real, use_banded)

    class(sysmtx_t), intent(inout)   :: this
    type(ext_complex_t), intent(out) :: det
    logical, optional, intent(in)    :: use_real
    logical, optional, intent(in)    :: use_banded

    logical          :: use_real_
    logical          :: use_banded_
    type(ext_real_t) :: det_r

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    if(PRESENT(use_banded)) then
       use_banded_ = use_banded
    else
       use_banded_ = .FALSE.
    endif

    ! Calculate the sysmtx determinant

    if(use_real_) then
       if(use_banded_) then
          call determinant_banded_r_(this, det_r)
       else
          call determinant_slu_r_(this, det_r)
       endif
       det = ext_complex_t(det_r)
    else
       if(use_banded_) then
          call determinant_banded_c_(this, det)
       else
          call determinant_slu_c_(this, det)
       end if
    endif

    ! Finish
    
    return

  end subroutine determinant_

!****

  $define $DETERMINANT_SLU $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine determinant_slu_${INFIX}_ (sm, det)

    class(sysmtx_t), intent(inout)   :: sm
    type(ext_${TYPE}_t), intent(out) :: det

    integer   :: n_e
    integer   :: n_i
    $TYPE(WP) :: M(2*sm%n_e,2*sm%n_e)
    integer   :: ipiv(2*sm%n_e)
    integer   :: info
    integer   :: i

    ! Factorize the matrix

    call factorize_slu_${INFIX}_(sm, det)

    ! Set up the reduced 2x2-block matrix

    n_e = sm%n_e
    n_i = sm%n_i

    $if ($INFIX eq 'r')

    M(:n_i,:n_e) = REAL(sm%B_i)
    M(n_i+1:n_i+n_e,:n_e) = REAL(sm%E_l(:,:,1))
    M(n_i+n_e+1:,:n_e) = 0._WP

    M(:n_i,n_e+1:) = 0._WP
    M(n_i+1:n_i+n_e,n_e+1:) = REAL(sm%E_r(:,:,1))
    M(n_i+n_e+1:,n_e+1:) = REAL(sm%B_o)

    $else

    M(:n_i,:n_e) = sm%B_i
    M(n_i+1:n_i+n_e,:n_e) = sm%E_l(:,:,1)
    M(n_i+n_e+1:,:n_e) = 0._WP

    M(:n_i,n_e+1:) = 0._WP
    M(n_i+1:n_i+n_e,n_e+1:) = sm%E_r(:,:,1)
    M(n_i+n_e+1:,n_e+1:) = sm%B_o

    $endif

    ! Add in its contribution to the determinant

    call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
    $ASSERT(info >= 0, Negative return from XGETRF)

    $if ($INFIX eq 'r')
    det = product([ext_real_t(diagonal(M)),det,ext_real_t(sm%S_i),ext_real_t(sm%S),ext_real_t(sm%S_o)])
    $else
    det = product([ext_complex_t(diagonal(M)),det,sm%S_i,sm%S,sm%S_o])
    $endif

    do i = 1,2*n_e
       if(ipiv(i) /= i) det = -det
    end do

    ! Finish

    return

  end subroutine determinant_slu_${INFIX}_

  $endsub

  $DETERMINANT_SLU(r,real)
  $DETERMINANT_SLU(c,complex)

!****

  $define $DETERMINANT_BANDED $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine determinant_banded_${INFIX}_ (sm, det)

    class(sysmtx_t), intent(inout)   :: sm
    type(ext_${TYPE}_t), intent(out) :: det

     $TYPE(WP), allocatable :: A_b(:,:)
    integer                 :: n_l
    integer                 :: n_u
    integer, allocatable    :: ipiv(:)
    integer                 :: info
    integer                 :: j

    ! Pack the smatrix into banded form

    call pack_banded_${INFIX}_(sm, A_b)

    ! LU decompose it

    n_l = sm%n_e + sm%n_i - 1
    n_u = sm%n_e + sm%n_i - 1

    allocate(ipiv(SIZE(A_b, 2)))

    $block

    $if ($DOUBLE_PRECISION)
    $if ($INFIX eq 'r')
    $local $X D
    $else
    $local $X Z
    $endif
    $else
    $if ($INFIX eq 'r')
    $local $X S
    $else
    $local $X C
    $endif
    $endif

    call ${X}GBTRF(SIZE(A_b, 2), SIZE(A_b, 2), n_l, n_u, A_b, SIZE(A_b, 1), ipiv, info)
    $ASSERT(info == 0 .OR. info == SIZE(A_b,2),Non-zero return from XGBTRF)

    $endblock

    ! Calculate the determinant

    $if ($INFIX eq 'r')
    det = product([ext_real_t(A_b(n_l+n_u+1,:)),ext_real_t(sm%S_i),ext_real_t(sm%S),ext_real_t(sm%S_o)])
    $else
    det = product([ext_complex_t(A_b(n_l+n_u+1,:)),sm%S_i,sm%S,sm%S_o])
    $endif
 
    do j = 1,SIZE(A_b, 2)
       if(ipiv(j) /= j) det = -det
    enddo

    ! Finish

    return

  end subroutine determinant_banded_${INFIX}_

  $endsub

  $DETERMINANT_BANDED(r,real)
  $DETERMINANT_BANDED(c,complex)

!****

  subroutine null_vector_ (this, b, det, use_real, use_banded)

    class(sysmtx_t), intent(inout)   :: this
    complex(WP), intent(out)         :: b(:)
    type(ext_complex_t), intent(out) :: det
    logical, optional, intent(in)    :: use_real
    logical, optional, intent(in)    :: use_banded

    logical          :: use_real_
    logical          :: use_banded_
    real(WP)         :: b_r(SIZE(b))
    type(ext_real_t) :: det_r

    $CHECK_BOUNDS(SIZE(b),this%n_e*(this%n+1))

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    if(PRESENT(use_banded)) then
       use_banded_ = use_banded
    else
       use_banded_ = .FALSE.
    endif

    ! Calculate the null vector

    if(use_real_) then
       if(use_banded_) then
          call null_vector_banded_r_(this, b_r, det_r)
       else
          call null_vector_slu_r_(this, b_r, det_r)
       endif
       b = b_r
       det = ext_complex_t(det_r)
    else
       if(use_banded_) then
          call null_vector_banded_c_(this, b, det)
       else
          call null_vector_slu_c_(this, b, det)
       endif
    end if

    ! Finish

    return

  end subroutine null_vector_

 !****

  $define $NULL_VECTOR_SLU $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine null_vector_slu_${INFIX}_ (sm, b, det)

    class(sysmtx_t), intent(inout)   :: sm
    $TYPE(WP), intent(out)           :: b(:)
    type(ext_${TYPE}_t), intent(out) :: det

    $TYPE(WP), parameter :: ZERO = 0._WP
    $TYPE(WP), parameter :: ONE = 1._WP

    integer   :: n
    integer   :: n_e
    integer   :: n_i
    $TYPE(WP) :: M(2*sm%n_e,2*sm%n_e)
    integer   :: ipiv(2*sm%n_e)
    integer   :: info
    $TYPE(WP) :: D(2*sm%n_e)
    integer   :: i
    $TYPE(WP) :: b_bound(2*sm%n_e)
    integer   :: l
    integer   :: k
    integer   :: i_a
    integer   :: i_b
    integer   :: i_c

    $CHECK_BOUNDS(SIZE(b),sm%n_e*(sm%n+1))

    ! Factorize the matrix

    call factorize_slu_${INFIX}_(sm, det, recon_mtx=.TRUE.)

    ! Set up the reduced 2x2-block matrix

    n = sm%n
    n_e = sm%n_e
    n_i = sm%n_i

    $if ($INFIX eq 'r')

    M(:n_i,:n_e) = REAL(sm%B_i)
    M(n_i+1:n_i+n_e,:n_e) = REAL(sm%E_l(:,:,1))
    M(n_i+n_e+1:,:n_e) = 0._WP

    M(:n_i,n_e+1:) = 0._WP
    M(n_i+1:n_i+n_e,n_e+1:) = REAL(sm%E_r(:,:,1))
    M(n_i+n_e+1:,n_e+1:) = REAL(sm%B_o)

    $else

    M(:n_i,:n_e) = sm%B_i
    M(n_i+1:n_i+n_e,:n_e) = sm%E_l(:,:,1)
    M(n_i+n_e+1:,:n_e) = 0._WP

    M(:n_i,n_e+1:) = 0._WP
    M(n_i+1:n_i+n_e,n_e+1:) = sm%E_r(:,:,1)
    M(n_i+n_e+1:,n_e+1:) = sm%B_o

    $endif

    ! Add in its contribution to the determinant

    call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
    $ASSERT(info >= 0, Negative return from XGETRF)

    $if ($INFIX eq 'r')
    det = product([ext_real_t(diagonal(M)),det,ext_real_t(sm%S)])
    $else
    det = product([ext_complex_t(diagonal(M)),det,sm%S])
    $endif

    do i = 1,2*n_e
       if(ipiv(i) /= i) det = -det
    end do

    ! Locate the singular element on the diagonal

    D = diagonal(M)

    if(ANY(D(:2*n_e-1) == 0._WP)) then
       i = MINLOC(ABS(D(:2*n_e-1)), DIM=1)
       $WARN(Singular element not at outer point)
    else
       i = 2*n_e
    endif

    ! Calculate the solutions at the two boundaries

    b_bound(:i-1) = -M(:i-1,i)

    $block

    $if ($DOUBLE_PRECISION)
    $if ($INFIX eq 'r')
    $local $X D
    $else
    $local $X Z
    $endif
    $else
    $if ($INFIX eq 'r')
    $local $X S
    $else
    $local $X C
    $endif
    $endif

    call ${X}TRSM('L', 'U', 'N', 'N', i-1, 1, &
         ONE, M, 2*n_e, b_bound, 2*n_e)

    b_bound(i) = 1._WP
    b_bound(i+1:) = 0._WP

    ! Now construct the full solution vector

    b(:n_e) = b_bound(:n_e)
    b(n_e*n+1:) = b_bound(n_e+1:)

    l = 1

    do
       if(l >= n) exit
       l = 2*l
    end do

    recon_loop : do

       l = l/2

       if (l == 0) exit recon_loop

       !$OMP PARALLEL DO SCHEDULE (DYNAMIC) PRIVATE (i_a, i_b, i_c)
       backsub_loop : do k = 1, n-l, 2*l

          i_a = (k-1)*n_e + 1
          i_b = i_a + l*n_e
          i_c = MIN(i_b + l*n_e, n_e*n+1)

          b(i_b:i_b+n_e-1) = ZERO

          $if ($INFIX eq 'r')
          call ${X}GEMV('N', n_e, n_e, -ONE, REAL(sm%E_l(:,:,k+l)), n_e, b(i_a:i_a+n_e-1), 1, ZERO, b(i_b:i_b+n_e-1), 1)
          call ${X}GEMV('N', n_e, n_e, -ONE, REAL(sm%E_r(:,:,k+l)), n_e, b(i_c:i_c+n_e-1), 1, ONE, b(i_b:i_b+n_e-1), 1)
          $else
          call ${X}GEMV('N', n_e, n_e, -ONE, sm%E_l(:,:,k+l), n_e, b(i_a:i_a+n_e-1), 1, ZERO, b(i_b:i_b+n_e-1), 1)
          call ${X}GEMV('N', n_e, n_e, -ONE, sm%E_r(:,:,k+l), n_e, b(i_c:i_c+n_e-1), 1, ONE, b(i_b:i_b+n_e-1), 1)
          $endif

       end do backsub_loop

    end do recon_loop

    $endblock

    ! Finish

    return

  end subroutine null_vector_slu_${INFIX}_

  $endsub

  $NULL_VECTOR_SLU(r,real)
  $NULL_VECTOR_SLU(c,complex)

!****

  $define $FACTORIZE_SLU $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine factorize_slu_${INFIX}_ (sm, det, recon_mtx)

    class(sysmtx_t), intent(inout)   :: sm
    type(ext_${TYPE}_t), intent(out) :: det
    logical, optional, intent(in)    :: recon_mtx

    $TYPE(WP), parameter :: ONE = 1._WP

    logical             :: recon_mtx_
    integer             :: n
    integer             :: n_e
    integer             :: l
    integer             :: k
    $TYPE(WP)           :: M_G(2*sm%n_e,sm%n_e)
    $TYPE(WP)           :: M_U(2*sm%n_e,sm%n_e)
    $TYPE(WP)           :: M_E(2*sm%n_e,sm%n_e)
    integer             :: ipiv(sm%n_e)
    integer             :: info
    integer             :: i
    type(ext_${TYPE}_t) :: block_det(sm%n)

    if(PRESENT(recon_mtx)) then
       recon_mtx_ = recon_mtx
    else
       recon_mtx_ = .FALSE.
    endif

    ! Factorize the sysmtx using the cyclic structured (SLU) algorithm
    ! by Wright (1994). The factorization is done in place, with
    ! E_l(:,:,1) and E_r(:,:,1) containing the final reduced blocks,
    ! and the other blocks of E_l and E_r (optionally) containing the
    ! U^-1 G and U^-1 E matrices needed to reconstruct solutions. The
    ! factorization determinant is also returned

    det = ext_${TYPE}_t(ONE)

    n = sm%n
    n_e = sm%n_e

    l = 1

    factor_loop : do

       if (l >= n) exit factor_loop

       ! Reduce pairs of blocks to single blocks

       !$OMP PARALLEL DO SCHEDULE (DYNAMIC) PRIVATE (M_G, M_U, M_E, ipiv, info, i)
       reduce_loop : do k = 1, n-l, 2*l

          ! Set up matrices (see expressions following eqn. 2.5 of
          ! Wright 1994)

          $if ($INFIX eq 'r')

          M_G(:n_e,:) = REAL(sm%E_l(:,:,k))
          M_G(n_e+1:,:) = 0._WP

          M_U(:n_e,:) = REAL(sm%E_r(:,:,k))
          M_U(n_e+1:,:) = REAL(sm%E_l(:,:,k+l))

          M_E(:n_e,:) = 0._WP
          M_E(n_e+1:,:) = REAL(sm%E_r(:,:,k+l))

          $else

          M_G(:n_e,:) = sm%E_l(:,:,k)
          M_G(n_e+1:,:) = 0._WP

          M_U(:n_e,:) = sm%E_r(:,:,k)
          M_U(n_e+1:,:) = sm%E_l(:,:,k+l)

          M_E(:n_e,:) = 0._WP
          M_E(n_e+1:,:) = sm%E_r(:,:,k+l)

          $endif

          ! Calculate the LU factorization of M_U, and use it to reduce
          ! M_E and M_G. The nasty fpx3 stuff is to ensure the correct
          ! LAPACK/BLAS routines are called (can't use generics, since
          ! we're then not allowed to pass array elements into
          ! assumed-size arrays; see, e.g., p. 268 of Metcalfe & Reid,
          ! "Fortran 90/95 Explained")

          call XGETRF(2*n_e, n_e, M_U, 2*n_e, ipiv, info)
          $ASSERT(info >= 0, Negative return from XGETRF)

          $block

          $if ($DOUBLE_PRECISION)
          $if ($INFIX eq 'r')
          $local $X D
          $else
          $local $X Z
          $endif
          $else
          $if ($INFIX eq 'r')
          $local $X S
          $else
          $local $X C
          $endif
          $endif

          call ${X}LASWP(n_e, M_E, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                        ONE, M_U(1,1), 2*n_e, M_E(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, -ONE, &
                        M_U(n_e+1,1), 2*n_e, M_E(1,1), 2*n_e, ONE, &
                        M_E(n_e+1,1), 2*n_e)

          if(recon_mtx_) then
             call ${X}TRSM('L', 'U', 'N', 'N', n_e, n_e, &
                           ONE, M_U(1,1), 2*n_e, M_E(1,1), 2*n_e)
          endif

          call ${X}LASWP(n_e, M_G, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                        ONE, M_U(1,1), 2*n_e, M_G(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, -ONE, &
                        M_U(n_e+1,1), 2*n_e, M_G(1,1), 2*n_e, ONE, &
                        M_G(n_e+1,1), 2*n_e)

          if(recon_mtx_) then
             call ${X}TRSM('L', 'U', 'N', 'N', n_e, n_e, &
                           ONE, M_U(1,1), 2*n_e, M_G(1,1), 2*n_e)
          endif

          $endblock

          ! Calculate the block determinant

          block_det(k) = product(ext_${TYPE}_t(diagonal(M_U)))

          do i = 1,n_e
             if(ipiv(i) /= i) block_det(k) = -block_det(k)
          end do

          ! Store results

          sm%E_l(:,:,k) = M_G(n_e+1:,:)
          sm%E_r(:,:,k) = M_E(n_e+1:,:)

          if(recon_mtx_) then
             sm%E_l(:,:,k+l) = M_G(:n_e,:)
             sm%E_r(:,:,k+l) = M_E(:n_e,:)
          endif

       end do reduce_loop

       ! Update the determinant

       det = product([block_det(:n-l:2*l),det])

       ! Loop around

       l = 2*l

    end do factor_loop

    ! Finish

    return

  end subroutine factorize_slu_${INFIX}_

  $endsub

  $FACTORIZE_SLU(r,real)
  $FACTORIZE_SLU(c,complex)

!****

  $define $PACK_BANDED $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine pack_banded_${INFIX}_ (sm, A_b)

    class(sysmtx_t), intent(in)         :: sm
    $TYPE(WP), allocatable, intent(out) :: A_b(:,:)

    integer :: n_l
    integer :: n_u
    integer :: k
    integer :: i_b
    integer :: j_b
    integer :: i
    integer :: j

    ! Pack the sysmtx into LAPACK's banded-matrix sparse format

    n_l = sm%n_e + sm%n_i - 1
    n_u = sm%n_e + sm%n_i - 1

    allocate(A_b(2*n_l+n_u+1,sm%n_e*(sm%n+1)))

    ! Inner boundary conditions
    
    A_b = 0._WP

    do j_b = 1, sm%n_e
       j = j_b
       do i_b = 1, sm%n_i
          i = i_b
          $if ($INFIX eq 'r')
          A_b(n_l+n_u+1+i-j,j) = REAL(sm%B_i(i_b,j_b))
          $else
          A_b(n_l+n_u+1+i-j,j) = sm%B_i(i_b,j_b)
          $endif
       end do
    end do

    ! Left equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = (k-1)*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + i_b + sm%n_i
             $if ($INFIX eq 'r')
             A_b(n_l+n_u+1+i-j,j) = REAL(sm%E_l(i_b,j_b,k))
             $else
             A_b(n_l+n_u+1+i-j,j) = sm%E_l(i_b,j_b,k)
             $endif
          end do
       end do
    end do

    ! Right equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = k*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + sm%n_i + i_b
             $if ($INFIX eq 'r')
             A_b(n_l+n_u+1+i-j,j) = REAL(sm%E_r(i_b,j_b,k))
             $else
             A_b(n_l+n_u+1+i-j,j) = sm%E_r(i_b,j_b,k)
             $endif
          end do
       end do
    end do

    ! Outer boundary conditions

    do j_b = 1, sm%n_e
       j = sm%n*sm%n_e + j_b
       do i_b = 1, sm%n_o
          i = sm%n*sm%n_e + sm%n_i + i_b
          $if ($INFIX eq 'r')
          A_b(n_l+n_u+1+i-j,j) = REAL(sm%B_o(i_b,j_b))
          $else
          A_b(n_l+n_u+1+i-j,j) = sm%B_o(i_b,j_b)
          $endif
       end do
    end do

    ! Finish

    return

  end subroutine pack_banded_${INFIX}_

  $endsub

  $PACK_BANDED(r,real)
  $PACK_BANDED(c,complex)

!****

  $define $NULL_VECTOR_BANDED $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine null_vector_banded_${INFIX}_ (sm, b, det)
  
    class(sysmtx_t), intent(inout)   :: sm
    $TYPE(WP), intent(out)           :: b(:)
    type(ext_${TYPE}_t), intent(out) :: det

    integer                :: m
    $TYPE(WP), allocatable :: A_b(:,:)
    integer                :: n_l
    integer                :: n_u
    integer, allocatable   :: ipiv(:)
    integer                :: info
    integer                :: i
    $TYPE(WP), allocatable :: A_r(:,:)
    $TYPE(WP), allocatable :: Mb(:,:)
    integer                :: j
    integer                :: n_lu

    $CHECK_BOUNDS(SIZE(b),sm%n_e*(sm%n+1))
    
    ! Pack the smatrix into banded form

    call pack_banded_${INFIX}_(sm, A_b)

    ! LU decompose it

    n_l = sm%n_e + sm%n_i - 1
    n_u = sm%n_e + sm%n_i - 1

    m = SIZE(A_b, 2)

    allocate(ipiv(m))

    $block

    $if ($DOUBLE_PRECISION)
    $if ($INFIX eq 'r')
    $local $X D
    $else
    $local $X Z
    $endif
    $else
    $if ($INFIX eq 'r')
    $local $X S
    $else
    $local $X C
    $endif
    $endif

    call ${X}GBTRF(m, m, n_l, n_u, A_b, SIZE(A_b, 1), ipiv, info)
    $ASSERT(info == 0 .OR. info == m,Non-zero return from XGBTRF)

    $endblock

    ! Calculate the determinant

    $if ($INFIX eq 'r')
    det = product([ext_real_t(A_b(n_l+n_u+1,:)),ext_real_t(sm%S)])
    $else
    det = product([ext_complex_t(A_b(n_l+n_u+1,:)),sm%S])
    $endif
 
    do j = 1, m
       if(ipiv(j) /= j) det = -det
    enddo

    ! Locate the singular element on the diagonal

    associate (D => A_b(n_l+n_u+1,:))

      if(ANY(D(:m-1) == 0._WP)) then
         i = MINLOC(ABS(D(:m-1)), DIM=1)
         $WARN(Singular element not at outer point)
      else
         i = m
      endif

    end associate

    ! Backsubstitute to solve the banded linear system A_b b = 0

    allocate(A_r(2*n_l+n_u+1,i-1))
    allocate(Mb(i-1,1))

    deallocate(ipiv)
    allocate(ipiv(i-1))

    if(i > 1) then 

       ! Set up the reduced LU system

       A_r(:n_l+n_u+1,:) = A_b(:n_l+n_u+1,:i-1)
       A_r(n_l+n_u+2:,:) = 0._WP

       ! The following line seems to cause out-of-memory errors when
       ! compiled with gfortran 4.8.0 on MVAPICH systems. Very puzzling!
       !
       ! ipiv = [(j,j=1,i-1)]

       do j = 1,i-1
          ipiv(j) = j
       enddo

       ! Solve for the 1:i-1 components of b

       n_lu = MIN(n_l+n_u, i-1)

       Mb(:i-n_lu-1,1) = 0._WP
       Mb(i-n_lu:,1) = -A_b(n_l+n_u+1-n_lu:n_l+n_u,i)

       $block

       $if ($DOUBLE_PRECISION)
       $if ($INFIX eq 'r')
       $local $X D
       $else
       $local $X Z
       $endif
       $else
       $if ($INFIX eq 'r')
       $local $X S
       $else
       $local $X C
       $endif
       $endif

       call ${X}GBTRS('N', SIZE(A_r, 2), n_l, n_u, 1, A_r, SIZE(A_r, 1), ipiv, Mb, SIZE(Mb, 1), info)
       $ASSERT(info == 0,Non-zero return from XGBTRS)

       $endblock

       b(:i-1) = Mb(:,1)

    end if
       
    ! Fill in the other parts of b
    
    b(i) = 1._WP
    b(i+1:) = 0._WP

    ! Finish

    return

  end subroutine null_vector_banded_${INFIX}_

  $endsub

  $NULL_VECTOR_BANDED(r,real)
  $NULL_VECTOR_BANDED(c,complex)

!****

  subroutine sol_vector_ (this, i, b)

    class(sysmtx_t), intent(inout)   :: this
    integer, intent(in)              :: i
    complex(WP), intent(out)         :: b(:)

    $CHECK_BOUNDS(SIZE(b),this%n_e*(this%n+1))

    ! Calculate the solution vector

    call sol_vector_banded_c_(this, i, b)

    ! Finish

    return

  end subroutine sol_vector_

!****

  $define $SOL_VECTOR_BANDED $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine sol_vector_banded_${INFIX}_ (sm, i, b)
  
    class(sysmtx_t), intent(inout)   :: sm
    integer, intent(in)              :: i
    $TYPE(WP), intent(out)           :: b(:)

    $TYPE(WP), allocatable :: A_b(:,:)
    integer                :: n_l
    integer                :: n_u
    integer                :: m
    integer, allocatable   :: ipiv(:)
    integer                :: info

    $CHECK_BOUNDS(SIZE(b),sm%n_e*(sm%n+1))
    
    ! Pack the smatrix into banded form

    call pack_banded_${INFIX}_(sm, A_b)

    ! Set up the RHS vector

    b = 0._WP
    b(i) = 1._WP

    ! Solve thge linear system

    n_l = sm%n_e + sm%n_i - 1
    n_u = sm%n_e + sm%n_i - 1

    m = SIZE(A_b, 2)

    allocate(ipiv(m))

    $block

    $if ($DOUBLE_PRECISION)
    $if ($INFIX eq 'r')
    $local $X D
    $else
    $local $X Z
    $endif
    $else
    $if ($INFIX eq 'r')
    $local $X S
    $else
    $local $X C
    $endif
    $endif

    call ${X}GBSV(m, n_l, n_u, 1, A_b, SIZE(A_b, 1), ipiv, b, SIZE(b, 1), info)
    $ASSERT(info == 0 .OR. info == m,Non-zero return from XGBSV)

    $endblock

    ! Finish

    return

  end subroutine sol_vector_banded_${INFIX}_

  $endsub

  $SOL_VECTOR_BANDED(r,real)
  $SOL_VECTOR_BANDED(c,complex)

end module gyre_sysmtx
