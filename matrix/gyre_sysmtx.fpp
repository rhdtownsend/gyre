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
     integer, allocatable             :: k_part(:)  ! Block partitioning indices
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
     procedure         :: determinant_banded
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

    allocate(this%k_part(OMP_SIZE_MAX+1))

    call partition_tasks(n, 2, this%k_part)

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

  recursive function determinant_slu (this) result (det)

    class(sysmtx_t), intent(inout) :: this
    type(ext_complex_t)            :: det

    integer             :: n_part(OMP_SIZE_MAX)
    integer             :: t
    type(ext_complex_t) :: elim_det(OMP_SIZE_MAX)
    type(sysmtx_t)      :: sm
    integer             :: k_src
    integer             :: k_dest

    ! Decide on the strategy

    if(this%n > 2) then

       ! Reduce the matrix using the structured factorization (SLU)
       ! algorithm by Wright (1994)

       n_part = this%k_part(2:) - this%k_part(:OMP_SIZE_MAX)

       ! Apply Gaussian elimination to the partitions

       !$OMP PARALLEL DO
       elim_loop : do t = 1,OMP_SIZE_MAX
          call elim_partition(this%E_l(:,:,this%k_part(t):this%k_part(t+1)-1), &
                              this%E_r(:,:,this%k_part(t):this%k_part(t+1)-1), &
                              elim_det=elim_det(t))
       end do elim_loop

       ! Initialize the reduced sysmtx

       call sm%init(COUNT(n_part > 0), this%n_e, this%n_i, this%n_o)

       sm%B_i = this%B_i
       sm%B_o = this%B_o

       ! Copy blocks to the reduced sysmtx

       !$OMP PARALLEL DO PRIVATE(k_src, k_dest)
       xfer_loop : do t = 1, OMP_SIZE_MAX

          if(n_part(t) > 0) then

             ! Determine the source and destination indices

             k_src = this%k_part(t+1)-1
             k_dest = COUNT(n_part(:t) > 0)

             ! Copy the blocks

             sm%E_l(:,:,k_dest) = this%E_l(:,:,k_src)
             sm%E_r(:,:,k_dest) = this%E_r(:,:,k_src)
             sm%S(k_dest) = product(this%S(this%k_part(t):this%k_part(t+1)-1))

          endif

       end do xfer_loop

       ! Combine the elimination determinants with the determinant of
       ! the reduced matrix

       det = product([elim_det,sm%determinant()])

    else

       ! Calculate the determinant using the banded algorithm

       det = this%determinant_banded()

    endif

    ! Finish

    return

  end function determinant_slu

!****

  function determinant_banded (this) result (det)

    class(sysmtx_t), intent(inout) :: this
    type(ext_complex_t)            :: det

    complex(WP), allocatable :: A_b(:,:)
    integer, allocatable     :: ipiv(:)
    integer                  :: n_l
    integer                  :: n_u
    integer                  :: i

    ! Pack the sysmtx into banded format

    call pack_banded(this, A_b)

    n_l = this%n_e + this%n_i - 1
    n_u = this%n_e + this%n_i - 1

    ! LU decompose the banded matrix

    allocate(ipiv(SIZE(A_b, 2)))

    call lu_decompose(A_b, n_l, n_u, ipiv)

    ! Calculate the determinant as the product of the diagonals (with
    ! sign flips to account for permutations)

    det = product(ext_complex(A_b(n_l+n_u+1,:)))

    do i = 1,SIZE(A_b, 2)
       if(ipiv(i) /= i) det = -det
    end do

    ! Add in the block scales

    det = product([det,this%S])

    ! Finish

    return

  end function determinant_banded

!****

  subroutine elim_partition (A, C, E, elim_det)

    complex(WP), intent(inout)                 :: A(:,:,:)
    complex(WP), intent(inout)                 :: C(:,:,:)
    complex(WP), intent(out), optional         :: E(:,:,:)
    type(ext_complex_t), intent(out), optional :: elim_det
      
    integer     :: n_e
    integer     :: n
    complex(WP) :: M_G(2*SIZE(A, 1),SIZE(A, 1))
    complex(WP) :: M_U(2*SIZE(A, 1),SIZE(A, 1))
    complex(WP) :: M_E(2*SIZE(A, 1),SIZE(A, 1))
    integer     :: k
    integer     :: ipiv(SIZE(A, 1))
    integer     :: info
    integer     :: i

    $ASSERT(SIZE(A, 2) == SIZE(A, 1),Dimension mismatch)

    $ASSERT(SIZE(C, 1) == SIZE(A, 1),Dimension mismatch)
    $ASSERT(SIZE(C, 2) == SIZE(A, 1),Dimension mismatch)
    $ASSERT(SIZE(C, 3) == SIZE(A, 3),Dimension mismatch)

    if(PRESENT(E)) then
       $ASSERT(SIZE(E, 1) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(E, 2) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(E, 3) == SIZE(A, 3)-1,Dimension mismatch)
    endif

    ! Perform the Gaussian elimination steps described in Section 2 of
    ! Wright (1994). A is overwritten by G and C by U

    n_e = SIZE(A, 1)
    n = SIZE(A, 3)

    if(PRESENT(elim_det)) elim_det = ext_complex(1._WP)

    if(n > 0) then

       M_G(n_e+1:,:) = A(:,:,1)
       M_E(n_e+1:,:) = C(:,:,1)

       block_loop : do k = 1,n-1
       
          ! Set up matrices (see expressions following eqn. 2.5 of
          ! Wright 1994)

          M_G(:n_e,:) = M_G(n_e+1:,:)
          M_G(n_e+1:,:) = 0._WP

          M_U(:n_e,:) = M_E(n_e+1:,:)
          M_U(n_e+1:,:) = A(:,:,k+1)

          M_E(:n_e,:) = 0._WP
          M_E(n_e+1:,:) = C(:,:,k+1)

          ! Transform the block by LU factoring M_U, and then
          ! multiplying the block by L^-1

          call LA_GETRF(2*n_e, n_e, M_U, 2*n_e, ipiv, info)
          $ASSERT(info == 0, Non-zero return from LA_GETRF)

          ! This nasty fpx3 stuff is to ensure the correct LAPACK/BLAS
          ! routines are called (can't use generics, since we're then
          ! not allowed to pass array elements into assumed-size
          ! arrays; see, e.g., p. 268 of Metcalfe & Reid, "Fortran
          ! 90/95 Explained")

          $block

          $if($DOUBLE_PRECISION)
          $local $X Z
          $else
          $local $X C
          $endif

          call ${X}LASWP(n_e, M_G, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                     CMPLX(1._WP, KIND=WP), M_U(1,1), 2*n_e, M_G(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(-1._WP, KIND=WP), &
                      M_U(n_e+1,1), 2*n_e, M_G(1,1), 2*n_e, CMPLX(1._WP, KIND=WP), &
                      M_G(n_e+1,1), 2*n_e)

          call ${X}LASWP(n_e, M_E, 2*n_e, 1, n_e, ipiv, 1)
          call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                     CMPLX(1._WP, KIND=WP), M_U(1,1), 2*n_e, M_E(1,1), 2*n_e)
          call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(-1._WP, KIND=WP), &
                      M_U(n_e+1,1), 2*n_e, M_E(1,1), 2*n_e, CMPLX(1._WP, KIND=WP), &
                      M_E(n_e+1,1), 2*n_e)

          $endblock

          ! Store results

          A(:,:,k) = M_G(:n_e,:)
          C(:,:,k) = M_U(:n_e,:)

          if(PRESENT(E)) E(:,:,k) = M_E(:n_e,:)

          ! Update the elimination determinant

          if(PRESENT(elim_det)) then

             elim_det = product([ext_complex(diagonal(M_U)),elim_det])

             do i = 1,n_e
                if(ipiv(i) /= i) elim_det = -elim_det
             end do
          
          endif

       end do block_loop

       ! Store final results

       A(:,:,n) = M_G(n_e+1:,:)
       C(:,:,n) = M_E(n_e+1:,:)

       ! Finish

       return

    endif

    ! Finish

    return

  end subroutine elim_partition

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
