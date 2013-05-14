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

  use gyre_ext_arith
  use gyre_linalg

  use f77_lapack

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
     type(ext_complex_t), allocatable :: S(:)       ! Block scales
     integer                          :: n          ! Number of equation blocks
     integer                          :: n_e        ! Number of equations per block
     integer                          :: n_i        ! Number of inner boundary conditions
     integer                          :: n_o        ! Number of outer boundary conditions
   contains
     private
     procedure, public :: init
     procedure, public :: set_inner_bound
     procedure, public :: set_outer_bound
     procedure, public :: set_block
     procedure, public :: determinant => determinant_slu
     procedure, public :: null_vector => null_vector_banded
  end type sysmtx_t

  ! Access specifiers

  private

  public :: sysmtx_t

  ! Procedures

contains

  subroutine init (this, n, n_e, n_i, n_o)

    class(sysmtx_t), intent(out) :: this
    integer, intent(in)          :: n
    integer, intent(in)          :: n_e
    integer, intent(in)          :: n_i
    integer, intent(in)          :: n_o

    ! Initialize the sysmtx

    allocate(this%E_l(n_e,n_e,n))
    allocate(this%E_r(n_e,n_e,n))

    allocate(this%B_i(n_i,n_e))
    allocate(this%B_o(n_o,n_e))

    allocate(this%S(n))

    this%n = n
    this%n_e = n_e
    this%n_i = n_i
    this%n_o = n_o

    ! Finish

    return

  end subroutine init

!****

  subroutine set_inner_bound (this, B_i)

    class(sysmtx_t), intent(inout)  :: this
    complex(WP), intent(in)         :: B_i(:,:)

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    ! Set the inner boundary conditions

    this%B_i = B_i

    ! Finish

    return

  end subroutine set_inner_bound

!****

  subroutine set_outer_bound (this, B_o)

    class(sysmtx_t), intent(inout)  :: this
    complex(WP), intent(in)         :: B_o(:,:)

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Set the outer boundary conditions

    this%B_o = B_o

    ! Finish

    return

  end subroutine set_outer_bound

!****

  subroutine set_block (this, k, E_l, E_r, S)

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

  end subroutine set_block

!****

  function determinant_slu (this) result (det)

    class(sysmtx_t), intent(in) :: this
    type(ext_complex_t)         :: det

    complex(WP)         :: A(this%n_e,this%n_e,this%n)
    complex(WP)         :: C(this%n_e,this%n_e,this%n)
    integer             :: n
    integer             :: n_e
    integer             :: k
    complex(WP)         :: P(2*this%n_e,this%n_e)
    complex(WP)         :: Q(2*this%n_e,this%n_e)
    complex(WP)         :: R(2*this%n_e,this%n_e)
    integer             :: ipiv(this%n_e)
    integer             :: info
    integer             :: i
    type(ext_complex_t) :: block_det(this%n)
    integer             :: n_red
    integer             :: n_i
    complex(WP)         :: M(2*this%n_e,2*this%n_e)
    integer             :: ipiv2(2*this%n_e)

    det = ext_complex(1._WP)

    ! Extract interior blocks from the sysmtx

    A = this%E_l
    C = this%E_r

    ! Repeatedly halve the number of these blocks by using the
    ! structured factorization (SLU) algorithm (Wright 1994) with a
    ! partition size of one or two

    n = this%n
    n_e = this%n_e

    factor_loop : do

       if (n==1) exit factor_loop

       ! Reduce pairs of blocks to single blocks

       !$OMP PARALLEL DO SCHEDULE (DYNAMIC) PRIVATE (P, Q, R, ipiv, info, i)
       reduce_loop : do k = 1, n-1, 2

          ! Set up matrices (see expressions following eqn. 2.5 of
          ! Wright 1994)

          P(:n_e,:) = C(:,:,k)
          P(n_e+1:,:) = A(:,:,k+1)

          Q(:n_e,:) = 0._WP
          Q(n_e+1:,:) = C(:,:,k+1)

          R(:n_e,:) = A(:,:,k)
          R(n_e+1:,:) = 0._WP

          ! Calculate the LU factorization of P, and use it to reduce
          ! Q and R. The nasty fpx3 stuff is to ensure the correct
          ! LAPACK/BLAS routines are called (can't use generics, since
          ! we're then not allowed to pass array elements into
          ! assumed-size arrays; see, e.g., p. 268 of Metcalfe & Reid,
          ! "Fortran 90/95 Explained")

          call LA_GETRF(2*n_e, n_e, P, 2*n_e, ipiv, info)
          $ASSERT(info == 0, Non-zero return from LA_GETRF)

          $block

          $if($DOUBLE_PRECISION)
          $local $X Z
          $else
          $local $X C
          $endif
          
          call ${X}LASWP(n_e, Q, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                        CMPLX(1._WP, KIND=WP), P(1,1), 2*n_e, Q(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(-1._WP, KIND=WP), &
                        P(n_e+1,1), 2*n_e, Q(1,1), 2*n_e, CMPLX(1._WP, KIND=WP), &
                        Q(n_e+1,1), 2*n_e)

          call ${X}LASWP(n_e, R, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                        CMPLX(1._WP, KIND=WP), P(1,1), 2*n_e, R(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(-1._WP, KIND=WP), &
                        P(n_e+1,1), 2*n_e, R(1,1), 2*n_e, CMPLX(1._WP, KIND=WP), &
                        R(n_e+1,1), 2*n_e)

          $endblock

          ! Calculate the block determinant

          block_det(k) = product(ext_complex(diagonal(P)))
          
          do i = 1,n_e
             if(ipiv(i) /= i) block_det(k) = -block_det(k)
          end do
             
          ! Store results

          C(:,:,k) = Q(n_e+1:,:)
          A(:,:,k) = R(n_e+1:,:)

       end do reduce_loop

       if(MOD(n, 2) /= 0) block_det(n) = ext_complex(1._WP)

       ! Update the determinant

       det = product([block_det(:n:2),det])

       ! Repack the reduced matrix

       n_red = (n+1)/2

       A(:,:,:n_red) = A(:,:,:n:2)
       C(:,:,:n_red) = C(:,:,:n:2)

       n = n_red

       ! Loop around

    end do factor_loop

    ! Process the final matrix, consisting of the boundary conditions
    ! plus a single internal block

    n_i = this%n_i

    M(:n_i,:n_e) = this%B_i
    M(n_i+1:n_i+n_e,:n_e) = A(:,:,1)
    M(n_i+n_e+1:,:n_e) = 0._WP

    M(:n_i,n_e+1:) = 0._WP
    M(n_i+1:n_i+n_e,n_e+1:) = C(:,:,1)
    M(n_i+n_e+1:,n_e+1:) = this%B_o

    call LA_GETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv2, info)
    $ASSERT(info == 0, Non-zero return from LA_GETRF)

    det = product([ext_complex(diagonal(M)),det,this%S])
    
    do i = 1,2*n_e
       if(ipiv2(i) /= i) det = -det
    end do

    ! Finish

    return

  end function determinant_slu

!****

  subroutine pack_banded (sm, A_b)

    class(sysmtx_t), intent(in)           :: sm
    complex(WP), allocatable, intent(out) :: A_b(:,:)

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
          A_b(n_l+n_u+1+i-j,j) = sm%B_i(i_b,j_b)
       end do
    end do

    ! Left equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = (k-1)*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + i_b + sm%n_i
             A_b(n_l+n_u+1+i-j,j) = sm%E_l(i_b,j_b,k)
          end do
       end do
    end do

    ! Right equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = k*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + sm%n_i + i_b
             A_b(n_l+n_u+1+i-j,j) = sm%E_r(i_b,j_b,k)
          end do
       end do
    end do

    ! Outer boundary conditions

    do j_b = 1, sm%n_e
       j = sm%n*sm%n_e + j_b
       do i_b = 1, sm%n_o
          i = sm%n*sm%n_e + sm%n_i + i_b
          A_b(n_l+n_u+1+i-j,j) = sm%B_o(i_b,j_b)
       end do
    end do

    ! Finish

    return

  end subroutine pack_banded

!****

  function null_vector_banded (this) result(b)

    class(sysmtx_t), intent(in) :: this
    complex(WP)                 :: b(this%n_e*(this%n+1))

    complex(WP), allocatable :: A_b(:,:)
    integer, allocatable     :: ipiv(:)
    integer                  :: n_l
    integer                  :: n_u
    integer                  :: i

    ! Pack the smatrix into banded form

    call pack_banded(this, A_b)

    ! LU decompose it

    allocate(ipiv(SIZE(A_b, 2)))

    n_l = this%n_e + this%n_i - 1
    n_u = this%n_e + this%n_i - 1

    call lu_decompose(A_b, n_l, n_u, ipiv)

    ! Locate the smallest diagonal element

    i = MINLOC(ABS(A_b(n_l+n_u+1,:)), DIM=1)

    if(SIZE(A_b, 2)-i > this%n_e) then
       $WARN(Smallest element not in final block)
    endif

    call lu_null_vector(A_b, n_l, n_u, i, b)

    ! Finish

    return

  end function null_vector_banded

end module gyre_sysmtx
