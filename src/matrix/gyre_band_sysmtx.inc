! Incfile  : gyre_band_sysmtx
! Purpose  : system matrix (banded storage, template)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_${T}_band_sysmtx

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

  ! Parameter definitions

  integer, parameter :: ALGO_LEN = 16

  ! Derived-type definitions

  type, extends (${T}_sysmtx_t) :: ${T}_sysmtx_t
     private
     $TYPE(WP), allocatable        :: A_b(:,:) ! Banded matrix
     integer, allocatable          :: ipiv(:)  ! Pivot indices
     type(${T}_ext_t)              :: S_i      ! Inner boundary scale
     type(${T}_ext_t)              :: S_o      ! Outer boundary scale
     type(${T}_ext_t), allocatable :: S(:)     ! Block scales
     integer                       :: n_ul     ! Number of sub-/super-diagonals
   contains
     private
     procedure, public :: set_B_i => set_B_i_
     procedure, public :: set_B_o => set_B_o_
     procedure, public :: set_block => set_block_
     procedure, public :: scale_rows => scale_rows_
     procedure, public :: determinant => determinant_
     procedure, public :: null_vector => null_vector_
     procedure, public :: sol_vector => sol_vector_
  end type ${T}_sysmtx_t

  ! Interfaces

  interface ${T}_sysmtx_t
     module procedure ${T}_sysmtx_t_
  end interface ${T}_sysmtx_t

  ! Access specifiers

  private

  public :: ${T}_sysmtx_t

  ! Procedures

contains

  function ${T}_band_sysmtx_t_ (n, n_e, n_i, n_o) result (sm)

    integer, intent(in)      :: n
    integer, intent(in)      :: n_e
    integer, intent(in)      :: n_i
    integer, intent(in)      :: n_o
    type(${T}_band_sysmtx_t) :: sm

    ! Construct the sysmtx_t

    ! Note that an additional n_ul rows are added to A_b to provide
    ! space for fill-in during factorization

    sm%n_ul = n_e + n_i - 1

    allocate(sm%A_b(3*sm%n_ul+1,n_e*(n+1)))
    allocate(sm%ipiv(n_e*(n+1)))

    allocate(sm%S(n))

    sm%n = n
    sm%n_e = n_e
    sm%n_i = n_i
    sm%n_o = n_o

    ! Finish

    return

  end function ${T}_sysmtx_t_

!****

  subroutine map_indices_ (this, i, j, i_b, j_b)

    class(${T}_band_sysmtx_t), intent(in) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(out) :: i_b
    integer, intent(out) :: j_b

    ! Map square matrix indices (i,j) to banded matrix indices
    ! (i_b,j_b); the latter are defined by LAPACK's banded storage
    ! format, with n_ul subdiagonals and 2*n_ul superdiagonals

    i_b = 2*this%n_ul + 1 + i - j
    j_b = j

    ! Finish

    return

  end subroutine map_indices_

!****

  subroutine set_B_i_ (this, B_i, S_i)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                    :: B_i(:,:)
    type(${T}_ext_t), intent(in)             :: S_i
    
    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    ! Set the inner boundary conditions

    ! Finish

    return

  end subroutine set_B_i_

!****

  subroutine set_B_o_ (this, B_o, S_o)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    $TYPE(WP), intent(in)                    :: B_o(:,:)
    type(${T}_ext_t), intent(in)             :: S_o

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Set the outer boundary conditions

    ! Finish

    return

  end subroutine set_B_o_

!****

  subroutine set_E_ (this, k, E_l, E_r, S)

    class(${T}_band_sysmtx_t), intent(inout) :: this
    integer, intent(in)                      :: k
    $TYPE(WP), intent(in)                    :: E_l(:,:)
    $TYPE(WP), intent(in)                    :: E_r(:,:)
    type(${T}_ext_t), intent(in)             :: S

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    $ASSERT(k >= 1,Invalid block index)
    $ASSERT(k <= this%n,Invalid block index)

    ! Set the block

    ! Finish

    return

  end subroutine set_block_

!****

  subroutine factorize_ (this)

    class(${T}_block_sysmtx_t), intent(inout) :: this

    $TYPE(WP), parameter :: ONE = 1._WP

    integer :: n
    integer :: info

    ! Factorize the sysmtx using LU decomposition

    n = SIZE(this%A_b, 1)
    m = SIZE(this%A_b, 2)

    call ${X}GBTRF(m, m, this%n_ul, this%n_ul, this%A_b, n, this%ipiv, info)
    $ASSERT(info == 0 .OR. info == m,Non-zero return from XGBTRF)

    ! Finish

    return

  end subroutine factorize_

!****

!   subroutine scale_rows_ (this)

!     class(${T}_sysmtx_t), intent(inout) :: this

!     real(WP) :: scale
!     integer  :: i
!     integer  :: k

!     ! Scale the rows of the sysmtx to have maximum absolute value of unity

!     do i = 1,this%n_i
!        scale = MAXVAL(ABS(this%B_i(i,:)))
!        this%B_i(i,:) = this%B_i(i,:)/scale
!        this%S_i = this%S_i*scale
!     end do

!     do k = 1, this%n
!        do i = 1, this%n_e
!           scale = MAX(MAXVAL(ABS(this%E_l(i,:,k))), MAXVAL(ABS(this%E_r(i,:,k))))
!           this%E_l(i,:,k) = this%E_l(i,:,k)/scale
!           this%E_r(i,:,k) = this%E_r(i,:,k)/scale
!           this%S(k) = this%S(k)*scale
!        end do
!     end do

!     do i = 1,this%n_o
!        scale = MAXVAL(ABS(this%B_o(i,:)))
!        this%B_o(i,:) = this%B_o(i,:)/scale
!        this%S_o = this%S_o*scale
!     end do

!     ! Finish

!     return

!   end subroutine scale_rows_

!****

  function det_ (this) result (det)

    class(${T}_band_sysmtx_t), intent(in) :: this
    type(${T}_ext_t)                      :: det

    integer :: j

    ! Evaluate the determinant

    det = product([${T}_ext_t(this%A_b(2*this%n_ul+1,:)),this%S_i,this%S,this%S_o])
 
    do j = 1, SIZE(this%ipiv)
       if(this%ipiv(j) /= j) det = -det
    enddo

    ! Finish

    return

  end function det_

!****

  function null_vec_ (this) result (v)
  
    class(${T}_band_sysmtx_t), intent(in) :: sm
    $TYPE(WP), intent(out)                :: v(this%n_e*(this%n+1))

    integer                :: i
    $TYPE(WP), allocatable :: A_r(:,:)
    $TYPE(WP), allocatable :: Mb(:,:)
    integer                :: j
    integer                :: n_lu

    ! Evaluate the null vector (working on the assumption that the
    ! final diagonal element of the factorized matrix is zero)

    associate (n => this%n, n_e => this%n_e, n_i => this%n_i)

      ! Locate the singular element on the diagonal

      sing_loop : do i = 1, n_e*(n+1)-1
         call this%map_indices_(i, i, i_b, j_b)
         if (this%A_b(i_b,j_b) == 0._WP)) exit sing_loop
      end do sing_loop

      if (i /= n_e*(n+1)) then
         $WARN(Singular element not at outer point)
      endif

      ! Calculate the solutions at the outer boundary

      v_bound(:i-1) = -M(:i-1,i)

      call ${X}TRSM('L', 'U', 'N', 'N', i-1, 1, &
           ONE, M, 2*n_e, v_bound, 2*n_e)

      v_bound(i) = 1._WP
      v_bound(i+1:) = 0._WP


       
      ! Reconstruct the full solution vector

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

         call ${X}GBTRS('N', SIZE(A_r, 2), n_l, n_u, 1, A_r, SIZE(A_r, 1), ipiv, Mb, SIZE(Mb, 1), info)
         $ASSERT(info == 0,Non-zero return from XGBTRS)

         b(:i-1) = Mb(:,1)

      end if
       
    ! Fill in the other parts of b
    
    b(i) = 1._WP
    b(i+1:) = 0._WP

    ! Finish

    return

  end subroutine null_vector_banded_

!****

  subroutine sol_vector_ (this, i, b)

    class(${T}_sysmtx_t), intent(inout) :: this
    integer, intent(in)                 :: i
    $TYPE(WP), intent(out)              :: b(:)

    $CHECK_BOUNDS(SIZE(b),this%n_e*(this%n+1))

    ! Calculate the solution vector

    call sol_vector_banded_(this, i, b)

    ! Finish

    return

  end subroutine sol_vector_

!****

  subroutine sol_vector_banded_ (sm, i, b)
  
    class(${T}_sysmtx_t), intent(inout) :: sm
    integer, intent(in)                 :: i
    $TYPE(WP), intent(out)              :: b(:)

    $TYPE(WP), allocatable :: A_b(:,:)
    integer                :: n_l
    integer                :: n_u
    integer                :: m
    integer, allocatable   :: ipiv(:)
    integer                :: info

    $CHECK_BOUNDS(SIZE(b),sm%n_e*(sm%n+1))
    
    ! Pack the smatrix into banded form

    call pack_banded_(sm, A_b)

    ! Set up the RHS vector

    b = 0._WP
    b(i) = 1._WP

    ! Solve thge linear system

    n_l = sm%n_e + sm%n_i - 1
    n_u = sm%n_e + sm%n_i - 1

    m = SIZE(A_b, 2)

    allocate(ipiv(m))

    call ${X}GBSV(m, n_l, n_u, 1, A_b, SIZE(A_b, 1), ipiv, b, SIZE(b, 1), info)
    $ASSERT(info == 0 .OR. info == m,Non-zero return from XGBSV)

    ! Finish

    return

  end subroutine sol_vector_banded_

end module gyre_${T}_sysmtx

$endsub

