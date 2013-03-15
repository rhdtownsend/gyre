! Module   : gyre_linalg
! Purpose  : linear algebra
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

module gyre_linalg

  ! Uses

  use core_kinds
  use core_order

  use f77_lapack
  use f95_lapack

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface gaussian_elim
     module procedure gaussian_elim_r
     module procedure gaussian_elim_c
  end interface gaussian_elim

  interface gaussian_mod
     module procedure gaussian_mod_r
     module procedure gaussian_mod_c
  end interface gaussian_mod

  interface lu_decompose
     module procedure lu_decompose_r
     module procedure lu_decompose_c
     module procedure lu_decompose_banded_r
     module procedure lu_decompose_banded_c
  end interface lu_decompose

  interface lu_backsub
     module procedure lu_backsub_r
     module procedure lu_backsub_c
     module procedure lu_backsub_banded_r
     module procedure lu_backsub_banded_c
  end interface lu_backsub

  interface lu_null_vector
     module procedure lu_null_vector_banded_r
     module procedure lu_null_vector_banded_c
  end interface lu_null_vector

  interface eigen_decompose
     module procedure eigen_decompose_r
     module procedure eigen_decompose_c
  end interface eigen_decompose

  interface linear_solve
     module procedure linear_solve_vec_r
     module procedure linear_solve_vec_c
     module procedure linear_solve_mat_r
     module procedure linear_solve_mat_c
  end interface linear_solve

  interface determinant
     module procedure determinant_r
     module procedure determinant_c
  end interface determinant

  interface commutator
     module procedure commutator_r
     module procedure commutator_c
  end interface commutator

  interface measure_bandwidth
     module procedure measure_bandwidth_full_r
     module procedure measure_bandwidth_full_c
  end interface measure_bandwidth

  interface diagonal
     module procedure diagonal_r
     module procedure diagonal_c
  end interface diagonal

  interface diagonal_matrix
     module procedure diagonal_matrix_r
     module procedure diagonal_matrix_c
  end interface diagonal_matrix

  interface matrix_exp
     module procedure matrix_exp_r
     module procedure matrix_exp_c
     module procedure matrix_exp_eigen
  end interface matrix_exp

  interface outer_product
     module procedure outer_product_r
     module procedure outer_product_c
     module procedure outer_product_l
  end interface outer_product

  interface partition
     module procedure partition_r
     module procedure partition_c
  end interface partition

  ! Access specifiers

  private

  public :: gaussian_elim
  public :: gaussian_mod
  public :: lu_decompose
  public :: lu_backsub
  public :: lu_null_vector
  public :: eigen_decompose
  public :: linear_solve
  public :: determinant
  public :: commutator
  public :: measure_bandwidth
  public :: diagonal
  public :: diagonal_matrix
  public :: identity_matrix
  public :: permutation_matrix
  public :: matrix_exp
  public :: outer_product
  public :: partition

  ! Procedures

contains

  $define $GAUSSIAN_ELIM $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine gaussian_elim_$SUFFIX (A, ipiv, scale_rows)
      
    $TYPE(WP), intent(inout)      :: A(:,:)
    integer, intent(out)          :: ipiv(:)
    logical, intent(in), optional :: scale_rows

    real(WP)  :: s(SIZE(A, 1))
    integer   :: k
    $TYPE(WP) :: A_tmp(SIZE(A, 2))
    real(WP)  :: s_tmp

    $ASSERT(SIZE(ipiv) == MIN(SIZE(A, 1), SIZE(A, 2)),Dimension mismatch)

    ! Perform Gaussian elimination on (generally, non-square) A, using
    ! algorithm 3.4.1 of Golub & Van Loan

    ! Determine row scale factors

    if(PRESENT(scale_rows)) then

       if(scale_rows) then

          s = MAXVAL(ABS(A), DIM=2)
          $ASSERT(ALL(s /= 0._WP),Singular matrix)

          s = 1._WP/s

       else

          s = 1._WP

       endif

    else

       s = 1._WP

    endif

    ! Loop down the diagonal

    do k = 1,MIN(SIZE(A, 1), SIZE(A, 2))

       ! Find the pivot row
         
       ipiv(k) = (k - 1) + MAXLOC(s(k:)*ABS(A(k:,k)), DIM=1)

       ! If required, swap rows

       if(ipiv(k) /= k) then

          A_tmp = A(k,:)
          A(k,:) = A(ipiv(k),:)
          A(ipiv(k),:) = A_tmp

          s_tmp = s(k)
          s(k) = s(ipiv(k))
          s(ipiv(k)) = s_tmp

       endif

       ! Update the matrix

       A(k+1:,k) = A(k+1:,k)/A(k,k)

       A(k+1:,k+1:) = A(k+1:,k+1:) - outer_product(A(k+1:,k), A(k,k+1:))

    end do

    ! Finish

    return

  end subroutine gaussian_elim_$SUFFIX

  $endsub

  $GAUSSIAN_ELIM(r,real)
  $GAUSSIAN_ELIM(c,complex)

!****

  $define $GAUSSIAN_MOD $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine gaussian_mod_$SUFFIX (A, ipiv, B)
      
    $TYPE(WP), intent(in)    :: A(:,:)
    integer, intent(in)      :: ipiv(:)
    $TYPE(WP), intent(inout) :: B(:,:)

    integer   :: k
    $TYPE(WP) :: B_tmp(SIZE(B, 2))
    integer   :: j

    $ASSERT(SIZE(ipiv) == MIN(SIZE(A, 1), SIZE(A, 2)),Dimension mismatch)

    $ASSERT(SIZE(B, 1) == SIZE(A, 1),Dimension mismatch)

    ! Perform Gaussian modification on B, using
    ! algorithm 3.4.1 of Golub & Van Loan

    ! Apply the permutation

    do k = 1,MIN(SIZE(A, 1), SIZE(A,2))

       if(ipiv(k) /= k) then

          B_tmp = B(k,:)
          B(k,:) = B(ipiv(k),:)
          B(ipiv(k),:) = B_tmp
          
       endif

    end do

    ! Do the forward sub

    do j = 1,SIZE(B, 2)

       do k = 1,MIN(SIZE(A, 1), SIZE(A, 2))

          ! Update the matrix

          B(k+1:,j) = B(k+1:,j) - B(k,j)*A(k+1:,k)

       end do

    end do

    ! Finish

    return

  end subroutine gaussian_mod_$SUFFIX

  $endsub

  $GAUSSIAN_MOD(r,real)
  $GAUSSIAN_MOD(c,complex)

!****

  $define $LU_DECOMPOSE $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine lu_decompose_$SUFFIX (A, ipiv, allow_singular)

    $TYPE(WP), intent(inout)       :: A(:,:)
    integer, intent(out), optional :: ipiv(:)
    logical, intent(in), optional  :: allow_singular

    logical :: allow_singular_
    integer :: ipiv_(SIZE(A, 1))
    integer :: info

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)

    if(PRESENT(ipiv)) then
       $ASSERT(SIZE(ipiv) == SIZE(A, 1),Dimension mismatch)
    endif

    if(PRESENT(allow_singular)) then
       allow_singular_ = allow_singular
    else
       allow_singular_ = .FALSE.
    endif

    ! LU decompose A

    call LA_GETRF(SIZE(A, 1), SIZE(A, 2), A, SIZE(A, 1), ipiv_, info)
    if(allow_singular_) then
       $ASSERT(info >= 0,Negative return from LA_GETRF)
    else
       $ASSERT(info == 0,Non-zero return from LA_GETRF)
    endif

    if(PRESENT(ipiv)) then
       ipiv = ipiv_
    endif

    ! Finish

    return

  end subroutine lu_decompose_$SUFFIX

  $endsub

  $LU_DECOMPOSE(r,real)
  $LU_DECOMPOSE(c,complex)

!****

  $define $LU_DECOMPOSE_BANDED $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine lu_decompose_banded_$SUFFIX (A_b, n_l, n_u, ipiv)

    $TYPE(WP), intent(inout)       :: A_b(:,:)
    integer, intent(in)            :: n_l
    integer, intent(in)            :: n_u
    integer, intent(out), optional :: ipiv(:)

    integer :: ipiv_(SIZE(A_b, 2))
    integer :: info

    $ASSERT(SIZE(A_b, 1) == 2*n_l+n_u+1,Dimension mismatch)

    if(PRESENT(ipiv)) then
       $ASSERT(SIZE(ipiv) == SIZE(A_b, 2),Dimension mismatch)
    endif

    ! LU decompose the banded (assumed-square) A_b

    call LA_GBTRF(SIZE(A_b, 2), SIZE(A_b, 2), n_l, n_u, A_b, SIZE(A_b, 1), ipiv_, info)
    $ASSERT(info == 0 .OR. info == SIZE(A_b,2),Non-zero return from LA_GBTRF)

    if(PRESENT(ipiv)) then
       ipiv = ipiv_
    endif

    ! Finish

    return

  end subroutine lu_decompose_banded_$SUFFIX

  $endsub

  $LU_DECOMPOSE_BANDED(r,real)
  $LU_DECOMPOSE_BANDED(c,complex)

!****

  $define $LU_BACKSUB $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine lu_backsub_$SUFFIX (A, ipiv, B)

    $TYPE(WP), intent(in)    :: A(:,:)
    integer, intent(in)      :: ipiv(:)
    $TYPE(WP), intent(inout) :: B(:,:)

    $TYPE(WP) :: A_(SIZE(A, 1),SIZE(A, 2))
    integer   :: info

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)
    $ASSERT(SIZE(ipiv) == SIZE(A, 1),Dimension mismatch)

    $ASSERT(SIZE(B, 1) == SIZE(A, 2),Dimension mismatch)

    ! Backsubstitute to solve the linear system A X = B

    A_ = A

    call LA_GETRS('N', SIZE(A_, 1), SIZE(B, 2), A_, SIZE(A_, 1), ipiv, B, SIZE(B, 1), info)
    $ASSERT(info == 0,Non-zero return from LA_GETRS)

    ! Finish

    return

  end subroutine lu_backsub_$SUFFIX

  $endsub

  $LU_BACKSUB(r, real)
  $LU_BACKSUB(c, complex)

!****

  $define $LU_BACKSUB_BANDED $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine lu_backsub_banded_$SUFFIX (A_b, n_l, n_u, ipiv, B)

    $TYPE(WP), intent(in)    :: A_b(:,:)
    integer, intent(in)      :: n_l
    integer, intent(in)      :: n_u
    integer, intent(in)      :: ipiv(:)
    $TYPE(WP), intent(inout) :: B(:,:)

    $TYPE(WP) :: A_b_(SIZE(A_b, 1),SIZE(A_b, 2))
    integer   :: info

    $ASSERT(SIZE(A_b, 1) == 2*n_l+n_u+1,Dimension mismatch)

    $ASSERT(SIZE(ipiv) == SIZE(A_b, 2),Dimension mismatch)

    $ASSERT(SIZE(B, 1) == SIZE(A_b, 2),Dimension mismatch)

    ! Backsubstitute to solve the banded linear system A_b B = B'

    A_b_ = A_b

    call LA_GBTRS('N', SIZE(A_b_, 2), n_l, n_u, SIZE(B, 2), A_b_, SIZE(A_b_, 1), ipiv, B, SIZE(B, 1), info)
    $ASSERT(info == 0,Non-zero return from LA_GBTRS)

    ! Finish

    return

  end subroutine lu_backsub_banded_$SUFFIX

  $endsub

  $LU_BACKSUB_BANDED(r, real)
  $LU_BACKSUB_BANDED(c, complex)

!****

  $define $LU_NULL_VECTOR_BANDED $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine lu_null_vector_banded_$SUFFIX (A_b, n_l, n_u, i, b)

    $TYPE(WP), intent(in)  :: A_b(:,:)
    integer, intent(in)    :: n_l
    integer, intent(in)    :: n_u
    integer, intent(in)    :: i
    $TYPE(WP), intent(out) :: b(:)

    $TYPE(WP) :: A_r(2*n_l+n_u+1,i-1)
    integer   :: j
    integer   :: ipiv(i-1)
    integer   :: n_lu
    integer   :: info

    $ASSERT(SIZE(A_b, 1) == 2*n_l+n_u+1,Dimension mismatch)

    $ASSERT(SIZE(b) == SIZE(A_b, 2),Dimension mismatch)

    ! Backsubstitute to solve the banded linear system A_b b = 0, for
    ! singular (rank = n-1) A (i.e., find the null vector of A_b). The
    ! singularity lies in the U(i,i) element of the (banded) LU
    ! decomposition of A; this element is not actually used, and
    ! assumed (without checking) to be zero

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

       b(:i-n_lu-1) = 0._WP
       b(i-n_lu:i-1) = -A_b(n_l+n_u+1-n_lu:n_l+n_u,i)

       call LA_GBTRS('N', SIZE(A_r, 2), n_l, n_u, 1, A_r, SIZE(A_r, 1), ipiv, b(:i-1), i-1, info)
       $ASSERT(info == 0,Non-zero return from LA_GBTRS)

    end if
       
    ! Fill in the other parts of b
    
    b(i) = 1._WP
    b(i+1:) = 0._WP

    ! Finish

    return

  end subroutine lu_null_vector_banded_$SUFFIX

  $endsub

  $LU_NULL_VECTOR_BANDED(r,real)
  $LU_NULL_VECTOR_BANDED(c,complex)

!****

  subroutine eigen_decompose_r (A, lambda, V_l, V_r, sort)

    real(WP), intent(in)               :: A(:,:)
    complex(WP), intent(out)           :: lambda(:)
    complex(WP), intent(out), optional :: V_l(:,:)
    complex(WP), intent(out), optional :: V_r(:,:)
    logical, intent(in), optional      :: sort

    logical     :: sort_
    real(WP)    :: A_(SIZE(A, 1),SIZE(A, 2))
    real(WP)    :: lambda_re(SIZE(lambda))
    real(WP)    :: lambda_im(SIZE(lambda))
    real(WP)    :: V_l_r(SIZE(A, 1),SIZE(A, 2))
    real(WP)    :: V_r_r(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_l_c(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_r_c(SIZE(A, 1),SIZE(A, 2))
    integer     :: ind(SIZE(A, 1))

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)
    $ASSERT(SIZE(lambda) == SIZE(A, 1),Dimension mismatch)

    if(PRESENT(V_l)) then
       $ASSERT(SIZE(V_l, 1) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(V_l, 2) == SIZE(A, 2),Dimension mismatch)
    endif
    
    if(PRESENT(V_r)) then
       $ASSERT(SIZE(V_r, 1) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(V_r, 2) == SIZE(A, 2),Dimension mismatch)
    endif

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    A_ = A

    if(PRESENT(V_l) .OR. PRESENT(V_r)) then

       call LA_GEEV(A_, lambda_re, lambda_im, VL=V_l_r, VR=V_r_r)

       call reconstruct_eigenvector(lambda_re, lambda_im, V_l_r, V_l_c)
       call reconstruct_eigenvector(lambda_re, lambda_im, V_r_r, V_r_c)

       V_l_c = CONJG(TRANSPOSE(V_l_c))

       ! Commented out because the LAPACK left/right eigenvectors
       ! don't seem to be proper inverses

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       ! do i = 1,SIZE(A,1)
       !    V_l_c(i,:) = V_l_c(i,:)/SUM(V_l_c(i,:)*V_r_c(:,i))
       ! enddo

       V_r_c = linear_solve(V_l_c, CMPLX(identity_matrix(SIZE(A, 1)), 0._WP, WP))

       if(PRESENT(V_l)) V_l = V_l_c
       if(PRESENT(V_r)) V_r = V_r_c

    else

       call LA_GEEV(A_, lambda_re, lambda_im)

    endif

    lambda = CMPLX(lambda_re, lambda_im, KIND=WP)

    ! (Possibly) sort

    if(sort_) then
       ind = sort_indices(REAL(lambda))
       lambda = lambda(ind)
       if(PRESENT(V_l)) V_l = V_l(ind,:)
       if(PRESENT(V_r)) V_r = V_r(:,ind)
    endif

    ! Finish

    return

  contains

    subroutine reconstruct_eigenvector (lambda_re, lambda_im, V, V_c)

      real(WP), intent(in)     :: lambda_re(:)
      real(WP), intent(in)     :: lambda_im(:)
      real(WP), intent(in)     :: V(:,:)
      complex(WP), intent(out) :: V_c(:,:)

      integer :: i

      ! Reconstruct the complex eigenvector

      i = 1
 
      do

         if(lambda_im(i) == 0._WP) then

            V_c(:,i) = V(:,i)

            i = i + 1

         elseif(lambda_re(i) == lambda_re(i+1) .AND. &
                lambda_im(i) == -lambda_im(i+1)) then

            V_c(:,i) = CMPLX(V(:,i), V(:,i+1), KIND=WP)
            V_c(:,i+1) = CONJG(V_c(:,i))

            i = i + 2

         else

            $ABORT(Complex eigenvalue does not have conjugate pair)

         endif

         if(i > SIZE(V, 2)) exit

      end do

    ! Finish

      return

    end subroutine reconstruct_eigenvector

  end subroutine eigen_decompose_r

!****

  subroutine eigen_decompose_c (A, lambda, V_l, V_r, sort)

    complex(WP), intent(in)            :: A(:,:)
    complex(WP), intent(out)           :: lambda(:)
    complex(WP), intent(out), optional :: V_l(:,:)
    complex(WP), intent(out), optional :: V_r(:,:)
    logical, intent(in), optional      :: sort

    logical     :: sort_
    complex(WP) :: A_(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_l_c(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_r_c(SIZE(A, 1),SIZE(A, 2))
    integer     :: ind(SIZE(A, 1))

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)
    $ASSERT(SIZE(lambda) == SIZE(A, 1),Dimension mismatch)

    if(PRESENT(V_l)) then
       $ASSERT(SIZE(V_l, 1) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(V_l, 2) == SIZE(A, 2),Dimension mismatch)
    endif
    
    if(PRESENT(V_r)) then
       $ASSERT(SIZE(V_r, 1) == SIZE(A, 1),Dimension mismatch)
       $ASSERT(SIZE(V_r, 2) == SIZE(A, 2),Dimension mismatch)
    endif

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    A_ = A

    if(PRESENT(V_l) .OR. PRESENT(V_r)) then

       call LA_GEEV(A_, lambda, VL=V_l_c, VR=V_r_c)

       V_l_c = CONJG(TRANSPOSE(V_l_c))

       ! Commented out because the LAPACK left/right eigenvectors
       ! don't seem to be proper inverses

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       ! do i = 1,SIZE(A, 1)
       !    V_l_c(i,:) = V_l_c(i,:)/SUM(V_l_c(i,:)*V_r_c(:,i))
       ! enddo

       V_r_c = linear_solve(V_l_c, CMPLX(identity_matrix(SIZE(A, 1)), 0._WP, WP))

       if(PRESENT(V_r)) V_r = V_r_c
       if(PRESENT(V_l)) V_l = V_l_c
       
    else

       call LA_GEEV(A_, lambda)

    endif

    ! (Possibly) sort

    if(sort_) then
       ind = sort_indices(REAL(lambda))
       lambda = lambda(ind)
       if(PRESENT(V_l)) V_l = V_l(ind,:)
       if(PRESENT(V_r)) V_r = V_r(:,ind)
    endif

    ! Finish

    return

  end subroutine eigen_decompose_c

!****

  $define $LINEAR_SOLVE $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function linear_solve_vec_$SUFFIX (A, b) result (x)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: b(:)
    $TYPE(WP)             :: x(SIZE(A, 2))

    $TYPE(WP) :: A_(SIZE(A, 1),SIZE(A, 2))

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)
    $ASSERT(SIZE(b) == SIZE(A, 1))

    ! Solve the linear system A x = b

    A_ = A
    x = b

    call LA_GESV(A_, x)

    ! Finish

    return

  end function linear_solve_vec_$SUFFIX

  function linear_solve_mat_$SUFFIX (A, B) result (X)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: B(:,:)
    $TYPE(WP)             :: X(SIZE(A, 2), SIZE(B, 2))

    $TYPE(WP) :: A_(SIZE(A, 1),SIZE(A, 2))

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)
    $ASSERT(SIZE(B, 1) == SIZE(A, 1))

    ! Solve the linear system A X = B

    A_ = A
    X = B

    call LA_GESV(A_, X)

    ! Finish

    return

  end function linear_solve_mat_$SUFFIX

  $endsub

  $LINEAR_SOLVE(r,real)
  $LINEAR_SOLVE(c,complex)

!****

  $define $DETERMINANT $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function determinant_$SUFFIX (A) result (det_A)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP)             :: det_A

    $TYPE(WP) :: A_(SIZE(A, 1),SIZE(A, 2))
    integer   :: ipiv(SIZE(A, 1))
    integer   :: i

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)

    ! LU decompose A

    A_ = A

    call lu_decompose(A_, ipiv, allow_singular=.TRUE.)

    ! Calculate the determinant from the diagonal

    det_A = 1._WP

    do i = 1,SIZE(A_, 1)
       det_A = det_A*A_(i,i)
       if(ipiv(i) /= i) det_A = -det_A
    end do

    ! Finish

    return

  end function determinant_$SUFFIX

  $endsub

  $DETERMINANT(r,real)
  $DETERMINANT(c,complex)

!****

  $define $COMMUTATOR $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function commutator_$SUFFIX (A, B) result (C)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: B(:,:)
    $TYPE(WP)             :: C(SIZE(A, 1),SIZE(A, 2))
    
    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)

    $ASSERT(SIZE(B, 1) == SIZE(B, 2),Dimension mismatch)
    $ASSERT(SIZE(B, 1) == SIZE(A, 1),Dimension mismatch)

    ! Calculate the commutator [A,B]

    C = MATMUL(A, B) - MATMUL(B, A)

    ! Finish

    return

  end function commutator_$SUFFIX

  $endsub

  $COMMUTATOR(r,real)
  $COMMUTATOR(c,complex)

!****

  $define $MEASURE_BANDWIDTH_FULL $sub
  
  $local $SUFFIX $1
  $local $TYPE $2
  $local $ZERO $3

  subroutine measure_bandwidth_full_$SUFFIX (A, n_l, n_u)

    $TYPE(WP), intent(in) :: A(:,:)
    integer, intent(out)  :: n_l
    integer, intent(out)  :: n_u

    integer :: n
    integer :: i
    integer :: j
    integer :: i_min
    integer :: i_max
    integer :: k

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)

    ! Measure the lower bandwidth

    n = SIZE(A, 1)

    lower_loop : do k = -(n-1), -1

       i_min = MAX(1, 1-k)
       i_max = MIN(SIZE(A, 1), SIZE(A, 1)-k)

       do i = i_min, i_max
          j = i + k
          if(A(i,j) /= $ZERO) exit lower_loop
       end do

    end do lower_loop

    n_l = -k

    ! Measure the upper bandwidth

    upper_loop : do k = n-1, 1, -1

       i_min = MAX(1, 1-k)
       i_max = MIN(SIZE(A, 1), SIZE(A, 1)-k)

       do i = i_min, i_max
          j = i + k
          if(A(i,j) /= $ZERO) exit upper_loop
       end do

    end do upper_loop

    n_u = k

    ! Finish

    return

  end subroutine measure_bandwidth_full_$SUFFIX

  $endsub

  $MEASURE_BANDWIDTH_FULL(r,real,0._WP)
  $MEASURE_BANDWIDTH_FULL(c,complex,(0._WP,0._WP))

!****

  $define $DIAGONAL $sub

  $local $SUFFIX $1
  $local $TYPE $2

  pure function diagonal_$SUFFIX (A) result (D)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP)             :: d(MIN(SIZE(A, 1),SIZE(A, 2)))

    integer :: i

    ! Extract the diagonal elements of matrix A

    do i = 1,SIZE(d)
       d(i) = A(i,i)
    end do

    ! Finish

    return

  end function diagonal_$SUFFIX

  $endsub

  $DIAGONAL(r,real)
  $DIAGONAL(c,complex)

!****

  $define $DIAGONAL_MATRIX $sub

  $local $SUFFIX $1
  $local $TYPE $2

  pure function diagonal_matrix_$SUFFIX (x) result (D)

    $TYPE(WP), intent(in) :: x(:)
    $TYPE(WP)             :: D(SIZE(x), SIZE(x))

    integer :: i

    ! Set up the diagonal matrix with elements x

    D = 0._WP

    do i = 1,SIZE(x)
       D(i,i) = x(i)
    end do

    ! Finish

    return

  end function diagonal_matrix_$SUFFIX

  $endsub

  $DIAGONAL_MATRIX(r,real)
  $DIAGONAL_MATRIX(c,complex)

!****

  pure function identity_matrix (n) result (I)

    integer, intent(in) :: n
    real(WP)            :: I(n,n)

    integer :: j

    ! Set up the rank-n identity matrix

    I = 0._WP

    do j = 1,n
       I(j,j) = 1._WP
    end do

    ! Finish

    return

  end function identity_matrix

!****

  pure function permutation_matrix (ipiv) result (P)

    integer, intent(in) :: ipiv(:)
    real(WP)            :: P(SIZE(ipiv),SIZE(ipiv))

    real(WP) :: tmp(SIZE(ipiv))
    integer  :: i
    
    ! Set up a permutation matrix from the row-exchange vector ipiv

    P = identity_matrix(SIZE(ipiv))

    do i = 1,SIZE(ipiv)
       if(ipiv(i) /= i) then
          tmp = P(i,:)
          P(i,:) = P(ipiv(i),:)
          P(ipiv(i),:) = tmp
       end if
    end do

    ! Finish

    return

  end function permutation_matrix

!****

  $define $MATRIX_EXP $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function matrix_exp_$SUFFIX (A, t) result (exp_At)

    $TYPE(WP), intent(in) :: A(:,:)
    real(WP), intent(in)  :: t
    $TYPE(WP)             :: exp_At(SIZE(A, 1),SIZE(A, 2))

    integer, parameter :: IDEG = 6

    complex(WP) :: lambda(SIZE(A,1))
    complex(WP) :: V_l(SIZE(A,1),SIZE(A,2))
    complex(WP) :: V_r(SIZE(A,1),SIZE(A,2))

    $ASSERT(SIZE(A, 1) == SIZE(A, 2),Dimension mismatch)

    ! Calculate the matrix exponential exp(A*t)

    call eigen_decompose_$SUFFIX(A, lambda, V_l=V_l, V_r=V_r)

    $if($TYPE eq "real")
    exp_At = REAL(matrix_exp_eigen(lambda, V_l, V_r, t))
    $else
    exp_At = matrix_exp_eigen(lambda, V_l, V_r, t)
    $endif

    ! Finish

    return

  end function matrix_exp_$SUFFIX

  $endsub

  $MATRIX_EXP(r,real)
  $MATRIX_EXP(c,complex)

!****

  function matrix_exp_eigen (lambda, V_l, V_r, t) result (exp_At)

    complex(WP), intent(in) :: lambda(:)
    complex(WP), intent(in) :: V_l(:,:)
    complex(WP), intent(in) :: V_r(:,:)
    real(WP), intent(in)    :: t
    complex(WP)             :: exp_At(SIZE(lambda),SIZE(lambda))

    integer :: i

    $ASSERT(SIZE(V_l, 2) == SIZE(V_l, 1),Dimension mismatch)
    $ASSERT(SIZE(V_r, 2) == SIZE(V_r, 1),Dimension mismatch)

    $ASSERT(SIZE(V_l, 1) == SIZE(lambda),Dimension mismatch)
    $ASSERT(SIZE(V_r, 1) == SIZE(lambda),Dimension mismatch)

    ! Calculate the matrix exponential exp(At) from the
    ! eigendecomposition of A

    do i = 1,SIZE(lambda)
       exp_At(i,:) = EXP(lambda(i)*t)*V_l(i,:)
    end do

    exp_At = MATMUL(V_r, exp_At)

    ! Finish

    return

  end function matrix_exp_eigen

!****

  $define $OUTER_PRODUCT $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function outer_product_$SUFFIX (v, w) result (A)

    $TYPE(WP), intent(in) :: v(:)
    $TYPE(WP), intent(in) :: w(:)
    $TYPE(WP)             :: A(SIZE(v), SIZE(w))

    ! Calculate the outer product v w^T

    A = SPREAD(v, 2, SIZE(w))*SPREAD(w, 1, SIZE(v))

    ! Finish

    return

  end function outer_product_$SUFFIX

  $endsub

  $OUTER_PRODUCT(r,real)
  $OUTER_PRODUCT(c,complex)

!****

  function outer_product_l (v, w) result (A)

    logical, intent(in) :: v(:)
    logical, intent(in) :: w(:)
    logical             :: A(SIZE(v), SIZE(w))

    ! Calculate the outer product v w^T

    A = SPREAD(v, 2, SIZE(w)) .AND. SPREAD(w, 1, SIZE(v))

    ! Finish

    return

  end function outer_product_l

!****

  $define $PARTITION $sub

  $local $SUFFIX $1
  $local $TYPE $2

  subroutine partition_$SUFFIX (A, mask_r, mask_c, A_part)

    $TYPE(WP), intent(in)  :: A(:,:)
    logical, intent(in)    :: mask_r(:)
    logical, intent(in)    :: mask_c(:)
    $TYPE(WP), intent(out) :: A_part(:,:)

    $ASSERT(SIZE(mask_r) == SIZE(A, 1),Dimension mismatch)
    $ASSERT(SIZE(mask_c) == SIZE(A, 2),Dimension mismatch)

    $ASSERT(COUNT(mask_r) == SIZE(A_part, 1),Dimension mismatch)
    $ASSERT(COUNT(mask_c) == SIZE(A_part, 2),Dimension mismatch)

    ! Partition A according to the row and column masks

    A_part = RESHAPE(PACK(A, outer_product(mask_r, mask_c)), SHAPE(A_part))

    ! Finish

    return

  end subroutine partition_$SUFFIX

  $endsub

  $PARTITION(r,real)
  $PARTITION(c,complex)

end module gyre_linalg
