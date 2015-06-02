! Incfile  : gyre_block_sysmtx
! Purpose  : system matrix (block storage, template)
!
! Copyright 2013-2015 Rich Townsend
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
     type(${T}_ext_t)              :: S_i        ! Inner boundary scale
     type(${T}_ext_t)              :: S_o        ! Outer boundary scale
     type(${T}_ext_t), allocatable :: S(:)       ! Block scales
   contains
     private
     procedure, public :: set_B_i => set_B_i_
     procedure, public :: set_B_o => set_B_o_
     procedure, public :: set_E => set_E_
     procedure, public :: factorize => factorize_
     procedure, public :: scale_rows_
     procedure, public :: det => det_
     procedure, public :: null_vec => null_vec_
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

    allocate(sm%S(n))

    sm%n = n
    sm%n_e = n_e
    sm%n_i = n_i
    sm%n_o = n_o

    ! Finish

    return

  end function ${T}_block_sysmtx_t_

!****

  subroutine set_B_i_ (this, B, S)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                     :: B(:,:)
    type(${T}_ext_t), intent(in)              :: S
    
    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    ! Set the inner boundary conditions

    this%B_i = B
    this%S_i = S

    ! Finish

    return

  end subroutine set_B_i_

!****

  subroutine set_B_o_ (this, B, S)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                     :: B(:,:)
    type(${T}_ext_t), intent(in)              :: S

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    ! Set the outer boundary conditions

    this%B_o = B
    this%S_o = S

    ! Finish

    return

  end subroutine set_B_o_

!****

  subroutine set_E_ (this, k, E_l, E_r, S)

    class(${T}_block_sysmtx_t), intent(inout) :: this
    integer, intent(in)                       :: k
    $TYPE(WP), intent(in)                     :: E_l(:,:)
    $TYPE(WP), intent(in)                     :: E_r(:,:)
    type(${T}_ext_t), intent(in)              :: S

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

  end subroutine set_E_

!****

  subroutine factorize_ (this)

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

            ! Fold the factorization determinant into S

            this%S(k) = this%S(k)*product(${T}_ext_t(diagonal(M_U)))

            do i = 1,n_e
               if (ipiv(i) /= i) this%S(k) = -this%S(k)
            end do

         end do reduce_loop

         ! Loop around

         l = 2*l

      end do factor_loop

    end associate

    ! Finish

    return

  end subroutine factorize_

!****

  subroutine scale_rows_ (this)

    class(${T}_block_sysmtx_t), intent(inout) :: this

    real(WP) :: scale
    integer  :: i
    integer  :: k

    ! Scale the rows of the sysmtx to have maximum absolute value of
    ! unity

    do i = 1, this%n_i
       scale = MAXVAL(ABS(this%B_i(i,:)))
       this%B_i(i,:) = this%B_i(i,:)/scale
       this%S_i = this%S_i*scale
    end do

    !$OMP PARALLEL DO PRIVATE (i, scale)
    do k = 1, this%n
       do i = 1, this%n_e
          scale = MAX(MAXVAL(ABS(this%E_l(i,:,k))), MAXVAL(ABS(this%E_r(i,:,k))))
          this%E_l(i,:,k) = this%E_l(i,:,k)/scale
          this%E_r(i,:,k) = this%E_r(i,:,k)/scale
          this%S(k) = this%S(k)*scale
       end do
    end do

    do i = 1, this%n_o
       scale = MAXVAL(ABS(this%B_o(i,:)))
       this%B_o(i,:) = this%B_o(i,:)/scale
       this%S_o = this%S_o*scale
    end do

    ! Finish

    return

  end subroutine scale_rows_

!****

  function det_ (this) result (det)

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

      det = product([${T}_ext_t(diagonal(M)),this%S_i,this%S,this%S_o])

      do i = 1,2*n_e
         if(ipiv(i) /= i) det = -det
      end do

    end associate

    ! Finish

    return

  end function det_

!****

  function null_vec_ (this) result (v)

    class(${T}_block_sysmtx_t), intent(in) :: this
    $TYPE(WP)                              :: v(this%n_e*(this%n+1))

    logical, parameter   :: USE_SVD = .FALSE.
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
    integer   :: l
    integer   :: k
    integer   :: i_a
    integer   :: i_b
    integer   :: i_c

    ! Evaluate the null vector (working on the assumption that the
    ! final diagonal element of the factorized matrix is zero)

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

         v_bound = V_H(2*n_e,:)

      else

         ! Use LU decomposition

         call XGETRF(2*n_e, 2*n_e, M, 2*n_e, ipiv, info)
         $ASSERT(info >= 0, Negative return from XGETRF)

         ! Locate the singular element on the diagonal of the outer
         ! block

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

      ! Backsubstitute to reconstruct the full solution vector

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

  end function null_vec_

end module gyre_${T}_block_sysmtx

$endsub

