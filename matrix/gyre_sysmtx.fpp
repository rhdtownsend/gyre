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
     type(ext_complex_t), allocatable :: scale(:)   ! Block scales
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

    allocate(this%scale(n))

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

  subroutine set_block (this, k, E_l, E_r, scale)

    class(sysmtx_t), intent(inout)  :: this
    integer, intent(in)             :: k
    complex(WP), intent(in)         :: E_l(:,:)
    complex(WP), intent(in)         :: E_r(:,:)
    type(ext_complex_t), intent(in) :: scale

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    $ASSERT(k >= 1,Invalid block index)
    $ASSERT(k <= this%n,Invalid block index)

    ! Set the block

    this%E_l(:,:,k) = E_l
    this%E_r(:,:,k) = E_r

    this%scale(k) = scale

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
                              elim_det(t))
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
             sm%scale(k_dest) = product(this%scale(this%k_part(t):this%k_part(t+1)-1))

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

    det = product([det,this%scale])

    ! Finish

    return

  end function determinant_banded

!****

  subroutine elim_partition (A, C, elim_det)

    complex(WP), intent(inout)       :: A(:,:,:)
    complex(WP), intent(inout)       :: C(:,:,:)
    type(ext_complex_t), intent(out) :: elim_det
      
    integer     :: n_e
    integer     :: n
    complex(WP) :: G(SIZE(A, 1),SIZE(A, 1),SIZE(A, 3))
    integer     :: k
    complex(WP) :: V(2*SIZE(A, 1),SIZE(A, 1))
    integer     :: i
    complex(WP) :: W(2*SIZE(A, 1),SIZE(A, 1))
    integer     :: ipiv(SIZE(A, 1))

    $ASSERT(SIZE(A, 2) == SIZE(A, 1),Dimension mismatch)

    $ASSERT(SIZE(C, 1) == SIZE(A, 1),Dimension mismatch)
    $ASSERT(SIZE(C, 2) == SIZE(A, 1),Dimension mismatch)

    ! Perform the Gaussian elimination steps described in Section 2 of
    ! Wright (1994)

    n_e = SIZE(A, 1)
    n = SIZE(A, 3)

    elim_det = ext_complex(1._WP)

    if(n > 0) then

       G(:,:,1) = A(:,:,1)

       block_loop : do k = 1, n-1
       
          ! Set up the current double-height block and Gaussian eliminate it

          V(:n_e,:) = C(:,:,k)
          V(n_e+1:,:) = A(:,:,k+1)

          call gaussian_elim(V, ipiv)

          ! Update the elimination determinant

          elim_det = product([ext_complex(diagonal(V)),elim_det])

          do i = 1,n_e
             if(ipiv(i) /= i) elim_det = -elim_det
          end do
          
          ! Update the appropriate blocks

          associate(R => A, E => C)

            ! R

            do i = 1,n_e
               R(i,:i-1,k) = 0._WP
               R(i,i:,k) = V(i,i:)
            end do
            
            ! C & E

            W(:n_e,:) = 0._WP
            W(n_e+1:,:) = C(:,:,k+1)

            call gaussian_mod(V, ipiv, W)

            E(:,:,k) = W(:n_e,:)
            C(:,:,k+1) = W(n_e+1:,:)

            ! G

            W(:n_e,:) = G(:,:,k)
            W(n_e+1:,:) = 0._WP

            call gaussian_mod(V, ipiv, W)

            G(:,:,k) = W(:n_e,:)
            G(:,:,k+1) = W(n_e+1:,:)

          end associate

       end do block_loop

       ! Extract the results
      
       A(:,:,n) = G(:,:,n)
       C(:,:,n) = C(:,:,n)

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

  subroutine pack_full (sm, A)

    class(sysmtx_t), intent(in) :: sm
    complex(WP), intent(out)    :: A(:,:)

    integer :: k
    integer :: i_b
    integer :: j_b
    integer :: i
    integer :: j

    $ASSERT(SIZE(A, 1) == sm%n_e*(sm%n+1),Dimension mismatch)
    $ASSERT(SIZE(A, 2) == sm%n_e*(sm%n+1),Dimension mismatch)

    ! Pack the sysmtx into a full matrix

    ! Inner boundary conditions
    
    A = 0._WP

    do j_b = 1, sm%n_e
       j = j_b
       do i_b = 1, sm%n_i
          i = i_b
          A(i,j) = sm%B_i(i_b,j_b)
       end do
    end do

    ! Left equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = (k-1)*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + i_b + sm%n_i
             A(i,j) = sm%E_l(i_b,j_b,k)
          end do
       end do
    end do

    ! Right equation blocks

    do k = 1, sm%n
       do j_b = 1, sm%n_e
          j = k*sm%n_e + j_b
          do i_b = 1, sm%n_e
             i = (k-1)*sm%n_e + sm%n_i + i_b
             A(i,j) = sm%E_r(i_b,j_b,k)
          end do
       end do
    end do

    ! Outer boundary conditions

    do j_b = 1, sm%n_e
       j = sm%n*sm%n_e + j_b
       do i_b = 1, sm%n_o
          i = sm%n*sm%n_e + sm%n_i + i_b
          A(i,j) = sm%B_o(i_b,j_b)
       end do
    end do

    ! Finish

    return

  end subroutine pack_full

!****

  function null_vector_inviter (this) result(b)

    class(sysmtx_t), intent(in) :: this
    complex(WP)                 :: b(this%n_e*(this%n+1))

    integer, parameter  :: MAX_ITER = 25
    real(WP), parameter :: EPS = 4._WP*EPSILON(0._WP)

    complex(WP), allocatable :: A_b(:,:)
    integer, allocatable     :: ipiv(:)
    integer                  :: n_l
    integer                  :: n_u
    real(WP)                 :: scale
    real(WP)                 :: tol
    integer                  :: i
    complex(WP)              :: x(this%n_e*(this%n+1),1)
    complex(WP)              :: x_prev(this%n_e*(this%n+1),1)
    integer                  :: j
    complex(WP)              :: mu

    ! **** NOTE: THIS ROUTINE IS CURRENTLY OUT-OF-ACTION; IT GIVES
    ! **** MUCH POORER RESULTS THAN null_vector_banded, AND IT'S NOT
    ! **** CLEAR WHY

    ! Pack the sysmtx into banded form

    call pack_banded(this, A_b)

    ! LU decompose it

    allocate(ipiv(SIZE(A_b, 2)))

    n_l = this%n_e + this%n_i - 1
    n_u = this%n_e + this%n_i - 1

    call lu_decompose(A_b, n_l, n_u, ipiv)

    ! Fix up small elements on the diagonal

    scale = MINVAL(ABS(A_b(n_l+n_u+1,:)), MASK=ABS(A_b(n_l+n_u+1,:)) /= 0._WP)

    where(ABS(A_b(n_l+n_u+1,:)) == 0._WP)
       A_b(n_l+n_u+1,:) = EPS*scale
    end where

    ! Use inverse iteration to converge on the null vector

    tol = EPS*SQRT(REAL(SIZE(x), WP))

    x = 1._WP

    iter_loop : do i = 1,MAX_ITER

       x_prev = x

       call lu_backsub(A_b, n_l, n_u, ipiv, x)

       ! Rescale x so that (i) the largest element is real, and (ii)
       ! the overall length is 1

       j = MAXLOC(ABS(x(:,1)), DIM=1)
       x = x/x(j,1)

       x = x/SQRT(DOT_PRODUCT(x(:,1), x(:,1)))

       ! Check for convergence

       mu = DOT_PRODUCT(x(:,1), x_prev(:,1))

       if(REAL(mu) < 0._WP) then
          x = -x
          mu = -mu
       endif

       if(ABS(mu - 1._WP) < tol) exit iter_loop

    end do iter_loop

    if(i > MAX_ITER) then
       write(ERROR_UNIT, *) 'x*x_proj = ', mu
       $WARN(Reached limit with inverse iterations)
    endif

    b = x(:,1)

    ! Finish

    ! return

  end function null_vector_inviter

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
