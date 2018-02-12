! Incfile  : gyre_block_sysmtx
! Purpose  : system matrix (block storage, template)
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

$if ($DOUBLE_PRECISION)
$if ($T eq 'r')
$local $X D
$else
$local $X Z
$endif
$else
$if ($T eq 'r')
$local $X S
$else
$local $X C
$endif
$endif

module gyre_${T}_block_sysmtx

  ! Uses

  use core_kinds
  use core_parallel
  use core_linalg

  use gyre_ext
  use gyre_linalg
  use gyre_sysmtx

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_sysmtx_t) :: ${T}_block_sysmtx_t
     private
     $TYPE(WP), allocatable        :: B_i(:,:)   ! Inner boundary conditions
     $TYPE(WP), allocatable        :: B_o(:,:)   ! Outer boundary conditions
     $TYPE(WP), allocatable        :: E_l(:,:,:) ! Left equation blocks
     $TYPE(WP), allocatable        :: E_r(:,:,:) ! Right equation blocks
     $TYPE(WP), allocatable        :: scl_i(:)   ! Inner boundary scales
     $TYPE(WP), allocatable        :: scl_o(:)   ! Outer boundary scales
     type(${T}_ext_t), allocatable :: scl(:)     ! Block scales
   contains
     private
     procedure, public :: set_B_i
     procedure, public :: set_B_o
     procedure, public :: set_E
     procedure, public :: factor
     procedure         :: scale_rows_
     procedure, public :: det
     procedure, public :: soln_vec_hom
     procedure, public :: soln_vec_inhom
     procedure         :: backsub_
     procedure, public :: resd_vec
  end type ${T}_block_sysmtx_t

  ! Interfaces

  interface ${T}_block_sysmtx_t
     module procedure ${T}_block_sysmtx_t_
  end interface ${T}_block_sysmtx_t

  ! Access specifiers

  private

  public :: ${T}_block_sysmtx_t

  ! Procedures

contains

  function ${T}_block_sysmtx_t_ (n, n_e, n_i, n_o) result (sm)

    integer, intent(in)       :: n
    integer, intent(in)       :: n_e
    integer, intent(in)       :: n_i
    integer, intent(in)       :: n_o
    type(${T}_block_sysmtx_t) :: sm

    $CHECK_BOUNDS(n_i+n_o,n_e)

    ! Construct the block_sysmtx_t

    allocate(sm%E_l(n_e,n_e,n))
    allocate(sm%E_r(n_e,n_e,n))

    allocate(sm%B_i(n_i,n_e))
    allocate(sm%B_o(n_o,n_e))

    allocate(sm%scl_i(n_i))
    allocate(sm%scl_o(n_o))

    allocate(sm%scl(n))

    sm%n = n
    sm%n_e = n_e
    sm%n_i = n_i
    sm%n_o = n_o

    ! Finish

    return

  end function ${T}_block_sysmtx_t_

  !****

  subroutine set_B_i (this, B, scl)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                     :: B(:,:)
    $TYPE(WP), intent(in)                     :: scl(:)
    
    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Set the inner boundary conditions

    this%B_i = B
    this%scl_i = scl

    ! Finish

    return

  end subroutine set_B_i

  !****

  subroutine set_B_o (this, B, scl)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                     :: B(:,:)
    $TYPE(WP), intent(in)                     :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Set the outer boundary conditions

    this%B_o = B
    this%scl_o = scl

    ! Finish

    return

  end subroutine set_B_o

  !****

  subroutine set_E (this, k, E_l, E_r, scl)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    integer, intent(in)                       :: k
    $TYPE(WP), intent(in)                     :: E_l(:,:)
    $TYPE(WP), intent(in)                     :: E_r(:,:)
    type(${T}_ext_t), intent(in)              :: scl

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    $ASSERT(k >= 1,Invalid block index)
    $ASSERT(k <= this%n,Invalid block index)

    ! Set the block

    this%E_l(:,:,k) = E_l
    this%E_r(:,:,k) = E_r

    this%scl(k) = scl

    ! Finish

    return

  end subroutine set_E

  !****

  subroutine factor (this)

    class(${T}_block_sysmtx_t), intent(inout) :: this

    $TYPE(WP), parameter :: ONE = 1._WP

    integer   :: l
    integer   :: k
    $TYPE(WP) :: M_G(2*this%n_e,this%n_e)
    $TYPE(WP) :: M_U(2*this%n_e,this%n_e)
    $TYPE(WP) :: M_E(2*this%n_e,this%n_e)
    integer   :: ipiv(this%n_e)
    integer   :: info
    integer   :: i

    ! Factorize the sysmtx using the cyclic structured (SLU) algorithm
    ! by [Wri1994]. The factorization is done in place: E_l(:,:,1) and
    ! E_r(:,:,1) contain the final reduced blocks, the other blocks of
    ! E_l and E_r contain the U^-1 G and U^-1 E matrices needed to
    ! reconstruct solutions, and S is scaled by the factorization
    ! determinants

    call this%scale_rows_()

    associate (n => this%n, n_e => this%n_e)

      ! Loop through factorization levels

      l = 1

      factor_loop : do

         if (l >= n) exit factor_loop

         ! Reduce pairs of blocks to single blocks

         !$OMP PARALLEL DO SCHEDULE (DYNAMIC) PRIVATE (M_G, M_U, M_E, ipiv, info, i)
         reduce_loop : do k = 1, n-l, 2*l

            ! Set up matrices (see expressions following eqn. 2.5 of
            ! Wright 1994)

            M_G(:n_e,:) = this%E_l(:,:,k)
            M_G(n_e+1:,:) = 0._WP

            M_U(:n_e,:) = this%E_r(:,:,k)
            M_U(n_e+1:,:) = this%E_l(:,:,k+l)

            M_E(:n_e,:) = 0._WP
            M_E(n_e+1:,:) = this%E_r(:,:,k+l)

            ! Calculate the LU factorization of M_U, and use it to reduce
            ! M_E and M_G

            call XGETRF(2*n_e, n_e, M_U, 2*n_e, ipiv, info)
            $ASSERT(info >= 0, Negative return from XGETRF)

            call ${X}LASWP(n_e, M_E, 2*n_e, 1, n_e, ipiv, 1)
            call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                          ONE, M_U(1,1), 2*n_e, M_E(1,1), 2*n_e)
            call ${X}GEMM('N', 'N', n_e, n_e, n_e, -ONE, &
                          M_U(n_e+1,1), 2*n_e, M_E(1,1), 2*n_e, ONE, &
                          M_E(n_e+1,1), 2*n_e)

            call ${X}TRSM('L', 'U', 'N', 'N', n_e, n_e, &
                          ONE, M_U(1,1), 2*n_e, M_E(1,1), 2*n_e)

            call ${X}LASWP(n_e, M_G, 2*n_e, 1, n_e, ipiv, 1)
            call ${X}TRSM('L', 'L', 'N', 'U', n_e, n_e, &
                 ONE, M_U(1,1), 2*n_e, M_G(1,1), 2*n_e)
            call ${X}GEMM('N', 'N', n_e, n_e, n_e, -ONE, &
                 M_U(n_e+1,1), 2*n_e, M_G(1,1), 2*n_e, ONE, &
                 M_G(n_e+1,1), 2*n_e)

            call ${X}TRSM('L', 'U', 'N', 'N', n_e, n_e, &
                          ONE, M_U(1,1), 2*n_e, M_G(1,1), 2*n_e)

            ! Store results

            this%E_l(:,:,k) = M_G(n_e+1:,:)
            this%E_r(:,:,k) = M_E(n_e+1:,:)

            this%E_l(:,:,k+l) = M_G(:n_e,:)
            this%E_r(:,:,k+l) = M_E(:n_e,:)

            ! Fold the factorization determinant into the scale

            this%scl(k) = this%scl(k)*PRODUCT(diagonal(M_U))

            do i = 1,n_e
               if (ipiv(i) /= i) this%scl(k) = -this%scl(k)
            end do

         end do reduce_loop

         ! Loop around

         l = 2*l

      end do factor_loop

    end associate

    ! Finish

    return

  end subroutine factor

  !****

  subroutine scale_rows_ (this)

    class(${T}_block_sysmtx_t), intent(inout) :: this

    real(WP) :: scl
    integer  :: i
    integer  :: k

    ! Scale the rows of the sysmtx to have maximum absolute value of
    ! unity

    do i = 1, this%n_i
       scl = MAXVAL(ABS(this%B_i(i,:)))
       this%B_i(i,:) = this%B_i(i,:)/scl
       this%scl_i(i) = this%scl_i(i)*scl
    end do

    !$OMP PARALLEL DO PRIVATE (i, scl)
    do k = 1, this%n
       do i = 1, this%n_e
          scl = MAX(MAXVAL(ABS(this%E_l(i,:,k))), MAXVAL(ABS(this%E_r(i,:,k))))
          this%E_l(i,:,k) = this%E_l(i,:,k)/scl
          this%E_r(i,:,k) = this%E_r(i,:,k)/scl
          this%scl(k) = this%scl(k)*scl
       end do
    end do

    do i = 1, this%n_o
       scl = MAXVAL(ABS(this%B_o(i,:)))
       this%B_o(i,:) = this%B_o(i,:)/scl
       this%scl_o(i) = this%scl_o(i)*scl
    end do

    ! Finish

    return

  end subroutine scale_rows_

  !****

  function det (this)

    class(${T}_block_sysmtx_t), intent(in) :: this
    type(${T}_ext_t)                       :: det

    $TYPE(WP) :: M(2*this%n_e,2*this%n_e)
    integer   :: ipiv(2*this%n_e)
    integer   :: info
    integer   :: i

    ! Evaluate the determinant

    associate (n_e => this%n_e, n_i => this%n_i)

      ! Set up the reduced 2x2-block matrix

      M(:n_i,:n_e) = this%B_i
      M(n_i+1:n_i+n_e,:n_e) = this%E_l(:,:,1)
      M(n_i+n_e+1:,:n_e) = 0._WP

      M(:n_i,n_e+1:) = 0._WP
      M(n_i+1:n_i+n_e,n_e+1:) = this%E_r(:,:,1)
      M(n_i+n_e+1:,n_e+1:) = this%B_o

      ! Factorize it

      call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
      $ASSERT(info >= 0, Negative return from XGETRF)

      ! Add together all the contributions to the determinant

      det = product([${T}_ext_t(diagonal(M)),${T}_ext_t(this%scl_i),this%scl,${T}_ext_t(this%scl_o)])

      do i = 1,2*n_e
         if(ipiv(i) /= i) det = -det
      end do

    end associate

    ! Finish

    return

  end function det

  !****

  function soln_vec_hom (this) result (v)

    class(${T}_block_sysmtx_t), intent(in) :: this
    $TYPE(WP)                              :: v(this%n_e*(this%n+1))

    logical, parameter   :: USE_SVD = .TRUE.
    $TYPE(WP), parameter :: ZERO = 0._WP
    $TYPE(WP), parameter :: ONE = 1._WP

    $TYPE(WP) :: M(2*this%n_e,2*this%n_e)
    real(WP)  :: sigma(2*this%n_e)
    $TYPE(WP) :: U(2*this%n_e,2*this%n_e)
    $TYPE(WP) :: V_H(2*this%n_e,2*this%n_e)
    integer   :: ipiv(2*this%n_e)
    integer   :: info
    integer   :: i_s
    real(WP)  :: M_s
    integer   :: i
    $TYPE(WP) :: v_bound(2*this%n_e)
 
    ! Evaluate the solution vector v of the homogeneous linear system
    ! S v = 0. It is assumed that the nullity nul(S) >= 1; if nul(S) =
    ! 0, then this routine will return the vector which minimizes the
    ! norm of S v. When nul(S) > 1, the routine returns the null
    ! vector associated with the smallest eigenvalue.

    associate (n => this%n, n_e => this%n_e, n_i => this%n_i)

      ! Set up the reduced 2x2-block matrix

      M(:n_i,:n_e) = this%B_i
      M(n_i+1:n_i+n_e,:n_e) = this%E_l(:,:,1)
      M(n_i+n_e+1:,:n_e) = 0._WP

      M(:n_i,n_e+1:) = 0._WP
      M(n_i+1:n_i+n_e,n_e+1:) = this%E_r(:,:,1)
      M(n_i+n_e+1:,n_e+1:) = this%B_o

      ! Solve for the solution at the two boundaries

      if (USE_SVD) then

         ! Use singular value decomposition

         call sing_decompose(M, sigma, U, V_H)

         $if ($T eq 'c')
         v_bound = CONJG(V_H(2*n_e,:))
         $else
         v_bound = V_H(2*n_e,:)
         $endif

      else

         ! Use LU decomposition

         call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
         $ASSERT(info >= 0, Negative return from XGETRF)

         ! Locate the smallest element on the diagonal of the outer
         ! block (this will be taken to be the singular element)

         i_s = 0
         M_s = HUGE(0._WP)

         sing_loop : do i = n_e+1, 2*n_e
            if (ABS(M(i,i)) < M_s) then
               M_s = ABS(M(i,i))
               i_s = i
            end if
         end do sing_loop

         v_bound(:i_s-1) = -M(:i_s-1,i_s)

         call ${X}TRSM('L', 'U', 'N', 'N', i_s-1, 1, &
                       ONE, M, 2*n_e, v_bound, 2*n_e)

         v_bound(i_s) = 1._WP
         v_bound(i_s+1:) = 0._WP

      endif

    end associate

    ! Backsubstitute to reconstruct the full solution vector

    call this%backsub_(v_bound, v)

    ! Finish

    return

  end function soln_vec_hom

  !****

  function soln_vec_inhom (this, w_i, w_o) result (v)

    class(${T}_block_sysmtx_t), intent(in) :: this
    $TYPE(WP), intent(in)                  :: w_i(:)
    $TYPE(WP), intent(in)                  :: w_o(:)
    $TYPE(WP)                              :: v(this%n_e*(this%n+1))

    $TYPE(WP) :: M(2*this%n_e,2*this%n_e)
    integer   :: ipiv(2*this%n_e)
    integer   :: info
    $TYPE(WP) :: B(2*this%n_e, 1)
    $TYPE(WP) :: v_bound(2*this%n_e)

    $CHECK_BOUNDS(SIZE(w_i),this%n_i)
    $CHECK_BOUNDS(SIZE(w_o),this%n_o)

    ! Evaluate the solution vector v of the inhomogeneous linear
    ! system S v = w. It is assumed that the right-hand side vector w
    ! has non-zero components in only the n_i first and n_o last rows
    ! (corresponding to the inner and outer boundary
    ! conditions). These components are supplied in w_i and w_o,
    ! respectively.
    
    associate (n => this%n, n_e => this%n_e, n_i => this%n_i)

      ! Set up the reduced 2x2-block matrix and rhs matrix

      M(:n_i,:n_e) = this%B_i
      M(n_i+1:n_i+n_e,:n_e) = this%E_l(:,:,1)
      M(n_i+n_e+1:,:n_e) = 0._WP

      M(:n_i,n_e+1:) = 0._WP
      M(n_i+1:n_i+n_e,n_e+1:) = this%E_r(:,:,1)
      M(n_i+n_e+1:,n_e+1:) = this%B_o

      B(:n_i,1) = w_i/this%scl_i
      B(n_i+1:n_i+n_e,1) = 0._WP
      B(n_i+n_e+1:,1) = w_o/this%scl_o

      ! Solve for the solution at the two boundaries using LU
      ! decomposition

      call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
      $ASSERT(info >= 0,Negative return from XGETRF)

      call XGETRS('N', 2*n_e, 1, M, 2*n_e, ipiv, B, 2*n_e, info)
      $ASSERT(info >= 0,Negative return from XGETRS)

      v_bound = B(:,1)

    end associate

    ! Backsubstitute to reconstruct the full solution vector

    call this%backsub_(v_bound, v)

    ! Finish

    return

  end function soln_vec_inhom

  !****

  subroutine backsub_ (this, v_bound, v)

    class(${T}_block_sysmtx_t), intent(in) :: this
    $TYPE(WP), intent(in)                  :: v_bound(:)
    $TYPE(WP), intent(out)                 :: v(:)

    $TYPE(WP), parameter :: ZERO = 0._WP
    $TYPE(WP), parameter :: ONE = 1._WP

    integer   :: l
    integer   :: k
    integer   :: i_a
    integer   :: i_b
    integer   :: i_c

    $CHECK_BOUNDS(SIZE(v_bound),2*this%n_e)
    $CHECK_BOUNDS(SIZE(v),this%n_e*(this%n+1))

    ! Given the solution vector v_bound at the boundaries,
    ! backsubstitute to reconstruct the full solution vector v

    associate (n => this%n, n_e => this%n_e, n_i => this%n_i)

      v(:n_e) = v_bound(:n_e)
      v(n_e*n+1:) = v_bound(n_e+1:)

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
            
            v(i_b:i_b+n_e-1) = ZERO

            call ${X}GEMV('N', n_e, n_e, -ONE, this%E_l(:,:,k+l), n_e, v(i_a:i_a+n_e-1), 1, ZERO, v(i_b:i_b+n_e-1), 1)
            call ${X}GEMV('N', n_e, n_e, -ONE, this%E_r(:,:,k+l), n_e, v(i_c:i_c+n_e-1), 1, ONE, v(i_b:i_b+n_e-1), 1)

         end do backsub_loop

      end do recon_loop

    end associate

    ! Finish

    return

  end subroutine backsub_

  !****

  function resd_vec (this, v) result (w)

    class(${T}_block_sysmtx_t), intent(in) :: this
    $TYPE(WP), intent(in)                  :: v(:)
    $TYPE(WP)                              :: w(this%n_e*(this%n+1))

    integer :: k
    integer :: i_a
    integer :: i_b

    $CHECK_BOUNDS(SIZE(v),this%n_e*(this%n+1))

    ! Evaluate the residual vector w = S v

    associate (n => this%n, n_e => this%n_e, n_i => this%n_i)

      w(:n_i) = MATMUL(this%B_i, v(:n_e))*this%scl_i

      !$OMP PARALLEL DO PRIVATE (i_a, i_b)
      multiply_loop : do k = 1, n

         i_a = (k-1)*n_e + 1
         i_b = i_a + n_i

         w(i_b:i_b+n_e-1) = MATMUL(this%E_l(:,:,k), v(i_a:i_a+n_e-1)) + MATMUL(this%E_r(:,:,k), v(i_a+n_e:i_a+2*n_e-1))

      end do multiply_loop

      w(n_e*n+n_i+1:) = MATMUL(this%B_o, v(n_e*n+1:))*this%scl_o

    end associate

    ! Finish

    return

  end function resd_vec
  
end module gyre_${T}_block_sysmtx

$endsub
