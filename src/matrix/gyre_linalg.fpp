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
  use core_linalg
  use core_order

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface eigen_decompose
     module procedure eigen_decompose_r_
     module procedure eigen_decompose_c_
  end interface eigen_decompose

  interface sing_decompose
     module procedure sing_decompose_r_
     module procedure sing_decompose_c_
  end interface sing_decompose

  interface linear_solve
     module procedure linear_solve_vec_r_
     module procedure linear_solve_vec_c_
     module procedure linear_solve_mat_r_
     module procedure linear_solve_mat_c_
  end interface linear_solve

  interface commutator
     module procedure commutator_r_
     module procedure commutator_c_
  end interface commutator

  interface measure_bandwidth
     module procedure measure_bandwidth_full_r_
     module procedure measure_bandwidth_full_c_
  end interface measure_bandwidth

  interface diagonal
     module procedure diagonal_r_
     module procedure diagonal_c_
  end interface diagonal

  interface diagonal_matrix
     module procedure diagonal_matrix_r_
     module procedure diagonal_matrix_c_
  end interface diagonal_matrix

  interface matrix_exp
     module procedure matrix_exp_r_
     module procedure matrix_exp_c_
     module procedure matrix_exp_eigen_
  end interface matrix_exp

  interface outer_product
     module procedure outer_product_r_
     module procedure outer_product_c_
     module procedure outer_product_l_
  end interface outer_product

  interface partition
     module procedure partition_r_
     module procedure partition_c_
  end interface partition

  ! Access specifiers

  private

  public :: eigen_decompose
  public :: sing_decompose
  public :: linear_solve
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

  subroutine eigen_decompose_r_ (A, lambda, V_l, V_r, sort)

    real(WP), intent(inout)       :: A(:,:)
    complex(WP), intent(out)      :: lambda(:)
    complex(WP), intent(out)      :: V_l(:,:)
    complex(WP), intent(out)      :: V_r(:,:)
    logical, optional, intent(in) :: sort

    logical     :: sort_
    integer     :: n
    real(WP)    :: lambda_re(SIZE(lambda))
    real(WP)    :: lambda_im(SIZE(lambda))
    real(WP)    :: V_l_r(SIZE(A, 1),SIZE(A, 2))
    real(WP)    :: V_r_r(SIZE(A, 1),SIZE(A, 2))
    real(WP)    :: work(4*SIZE(A, 1))
    integer     :: info
    integer     :: i
    complex(WP) :: norm
    integer     :: ind(SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(lambda),SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(A, 2))
    
    $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(A, 2))

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    n = SIZE(A, 1)

    if(n == 2) then

       ! Special case: n == 2

       call eigen_decompose_2_r_(A, lambda, V_l, V_r)

    else

       call XGEEV('V', 'V', n, A, n, lambda_re, lambda_im, &
                  V_l_r, n, V_r_r, n, work, SIZE(work), info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

       lambda = CMPLX(lambda_re, lambda_im, WP)

       call reconstruct_eigenvector_(lambda_re, lambda_im, V_l_r, V_l)
       call reconstruct_eigenvector_(lambda_re, lambda_im, V_r_r, V_r)

       V_l = CONJG(TRANSPOSE(V_l))

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       do i = 1, n
          norm = SUM(V_l(i,:)*V_r(:,i))
!          $ASSERT(ABS(norm) > n*EPSILON(0._WP),Near-defective matrix)
          V_l(i,:) = V_l(i,:)/norm
       enddo

    endif

    ! (Possibly) sort

    if(sort_) then

       ind = sort_indices(lambda_re)

       lambda = lambda(ind)
       
       V_l = V_l(ind,:)
       V_r = V_r(:,ind)

    endif

    ! Finish

    return

  contains

    subroutine reconstruct_eigenvector_ (lambda_re, lambda_im, V_r, V)

      real(WP), intent(in)     :: lambda_re(:)
      real(WP), intent(in)     :: lambda_im(:)
      real(WP), intent(in)     :: V_r(:,:)
      complex(WP), intent(out) :: V(:,:)

      integer :: i

      ! Reconstruct the complex eigenvector

      i = 1
 
      recon_loop : do

         if(lambda_im(i) == 0._WP) then

            V(:,i) = V_r(:,i)

            i = i + 1

         elseif(lambda_re(i) == lambda_re(i+1) .AND. &
                lambda_im(i) == -lambda_im(i+1)) then

            V(:,i) = CMPLX(V_r(:,i), V_r(:,i+1), KIND=WP)
            V(:,i+1) = CONJG(V(:,i))

            i = i + 2

         else

            $ABORT(Complex eigenvalue does not have conjugate pair)

         endif

         if(i > SIZE(V_r, 2)) exit recon_loop

      end do recon_loop

    ! Finish

      return

    end subroutine reconstruct_eigenvector_

  end subroutine eigen_decompose_r_

  !****

  subroutine eigen_decompose_c_ (A, lambda, V_l, V_r, sort)

    complex(WP), intent(inout)    :: A(:,:)
    complex(WP), intent(out)      :: lambda(:)
    complex(WP), intent(out)      :: V_l(:,:)
    complex(WP), intent(out)      :: V_r(:,:)
    logical, optional, intent(in) :: sort

    logical     :: sort_
    integer     :: n
    complex(WP) :: work(2*SIZE(A, 1))
    real(WP)    :: rwork(2*SIZE(A, 1))
    integer     :: info
    integer     :: i
    complex(WP) :: norm
    integer     :: ind(SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(lambda),SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(A, 2))
    
    $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(A, 2))

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    n = SIZE(A, 1)

    if(n == 2) then

       ! Special case: n == 2

       call eigen_decompose_2_c_(A, lambda, V_l, V_r)

    else

       call XGEEV('V', 'V', n, A, n, lambda, &
                   V_l, n, V_r, n, work, SIZE(work, 1), rwork, info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

       V_l = CONJG(TRANSPOSE(V_l))

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       do i = 1,n
          norm = SUM(V_l(i,:)*V_r(:,i))
!          $ASSERT(ABS(norm) > n*EPSILON(0._WP),Near-defective matrix)
          V_l(i,:) = V_l(i,:)/norm
       enddo

    end if

    ! (Possibly) sort

    if(sort_) then

       ind = sort_indices(REAL(lambda))

       lambda = lambda(ind)

       V_l = V_l(ind,:)
       V_r = V_r(:,ind)

    endif

    ! Finish

    return

  end subroutine eigen_decompose_c_

  !****

  $define $EIGEN_DECOMPOSE_2 $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine eigen_decompose_2_${INFIX}_ (A, lambda, V_l, V_r)

    $TYPE(WP), intent(in)              :: A(:,:)
    complex(WP), intent(out)           :: lambda(:)
    complex(WP), intent(out), optional :: V_l(:,:)
    complex(WP), intent(out), optional :: V_r(:,:)

    $TYPE(WP)   :: b
    $TYPE(WP)   :: c
    complex(WP) :: s
    complex(WP) :: q
    complex(WP) :: V_l_(2,2)
    complex(WP) :: V_r_(2,2)
    integer     :: i

    $CHECK_BOUNDS(SIZE(A, 1),2)
    $CHECK_BOUNDS(SIZE(A, 2),2)

    $CHECK_BOUNDS(SIZE(lambda),SIZE(A, 1))

    if(PRESENT(V_l)) then
       $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(A, 2))
    endif
    
    if(PRESENT(V_r)) then
       $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(A, 2))
    endif

    ! Perform the eigendecomposition of the 2x2 matrix A

    ! Solve the quadratic for the eigenvalues

    b = -A(1,1) - A(2,2)
    c = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    $if($TYPE eq 'real')
    s = SQRT(CMPLX(b**2 - 4._WP*c, KIND=KIND(0._WP)))
    q = -0.5_WP*(b + SIGN(1._WP, b)*s)
    $else
    s = SQRT(b**2 - 4._WP*c)
    q = -0.5_WP*(b + SIGN(1._WP, REAL(CONJG(b)*s))*s)
    $endif

    lambda(1) = q
    lambda(2) = c/q

    ! Calculate eigenvectors

    if(PRESENT(V_l) .OR. PRESENT(V_r)) then

       V_l_(1,1) = A(2,1)
       V_l_(1,2) = lambda(1) - A(1,1)

       V_l_(2,1) = A(2,1)
       V_l_(2,2) = lambda(2) - A(1,1)
    
       V_r_(1,1) = A(1,2)
       V_r_(2,1) = lambda(1) - A(1,1)
    
       V_r_(1,2) = A(1,2)
       V_r_(2,2) = lambda(2) - A(1,1)
       
       do i = 1,SIZE(A, 1)
          V_l_(i,:) = V_l_(i,:)/SUM(V_l_(i,:)*V_r_(:,i))
       enddo

       if(PRESENT(V_l)) then
          V_l = V_l_
       endif

       if(PRESENT(V_r)) then
          V_r = V_r_
       endif

    end if

    ! Finish

    return

  end subroutine eigen_decompose_2_${INFIX}_

  $endsub

  $EIGEN_DECOMPOSE_2(r,real)
  $EIGEN_DECOMPOSE_2(c,complex)

  !****

  subroutine sing_decompose_r_ (A, sigma, U, V_H)

    real(WP), intent(inout) :: A(:,:)
    real(WP), intent(out)   :: sigma(:)
    real(WP), intent(out)   :: U(:,:)
    real(WP), intent(out)   :: V_H(:,:)

    integer  :: n
    real(WP) :: work(5*SIZE(A, 1))
    integer  :: info

    $CHECK_BOUNDS(SIZE(sigma),MIN(SIZE(A, 1), SIZE(A,2)))

    $CHECK_BOUNDS(SIZE(U, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(U, 2),SIZE(A, 1))
    
    $CHECK_BOUNDS(SIZE(V_H, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(V_H, 2),SIZE(A, 2))

    ! Perform the singular-value decomposition of A

    n = SIZE(A, 1)

    call XGESVD('A', 'A', n, n, A, n, sigma, U, n, V_H, n, work, SIZE(work), info)
    $ASSERT(info == 0,Non-zero return from XGESVD)

    ! Finish

    return

  end subroutine sing_decompose_r_
    
  !****

  subroutine sing_decompose_c_ (A, sigma, U, V_H)

    complex(WP), intent(inout) :: A(:,:)
    real(WP), intent(out)      :: sigma(:)
    complex(WP), intent(out)   :: U(:,:)
    complex(WP), intent(out)   :: V_H(:,:)

    integer     :: n
    complex(WP) :: work(3*SIZE(A, 1))
    real(WP)    :: rwork(5*SIZE(A, 1))
    integer     :: info

    $CHECK_BOUNDS(SIZE(sigma),MIN(SIZE(A, 1), SIZE(A, 2)))

    $CHECK_BOUNDS(SIZE(U, 1),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(U, 2),SIZE(A, 1))
    
    $CHECK_BOUNDS(SIZE(V_H, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(V_H, 2),SIZE(A, 2))

    ! Perform the singular-value decomposition of A

    n = SIZE(A, 1)

    call XGESVD('A', 'A', n, n, A, n, sigma, U, n, V_H, n, work, SIZE(work), rwork, info)
    $ASSERT(info == 0,Non-zero return from XGESVD)

    ! Finish

    return

  end subroutine sing_decompose_c_
    
  !****

  $define $LINEAR_SOLVE $sub

  $local $INFIX $1
  $local $TYPE $2

  function linear_solve_vec_${INFIX}_ (A, b) result (x)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: b(:)
    $TYPE(WP)             :: x(SIZE(A, 2))

    $TYPE(WP) :: Mx(SIZE(x), 1)
    $TYPE(WP) :: Mb(SIZE(b), 1)

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(b),SIZE(A, 1))

    ! Solve the linear system A x = b

    Mb(:,1) = b

    Mx = linear_solve(A, Mb)

    x = Mx(:,1)

    ! Finish

    return

  end function linear_solve_vec_${INFIX}_

  function linear_solve_mat_${INFIX}_ (A, B) result (X)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: B(:,:)
    $TYPE(WP)             :: X(SIZE(A, 2), SIZE(B, 2))

    $TYPE(WP)    :: A_(SIZE(A, 1),SIZE(A, 2))
    $TYPE(WP)    :: B_(SIZE(B, 1),SIZE(B, 2))
    integer      :: ipiv(SIZE(A, 1))
    integer      :: info
    $TYPE(WP)    :: Af(SIZE(A, 1),SIZE(A, 2))
    real(WP)     :: R(SIZE(A, 1))
    real(WP)     :: C(SIZE(A, 1))
    real(WP)     :: rcond
    real(WP)     :: ferr(SIZE(B, 2))
    real(WP)     :: berr(SIZE(B, 2))
    $if($TYPE eq 'complex')
    $TYPE(WP)    :: work(2*SIZE(A, 1))
    real(WP)     :: xwork(2*SIZE(A, 1))
    $else
    $TYPE(WP)    :: work(4*SIZE(A, 1))
    integer      :: xwork(SIZE(A, 1))
    $endif
    character(1) :: equed

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(B, 1),SIZE(A, 1))

    ! Solve the linear system A X = B

    A_ = A
    B_ = B

    call XGESVX('E', 'N', SIZE(A_, 1), SIZE(B_, 2), A_, SIZE(A_, 1), AF, SIZE(AF, 1), ipiv, &
                equed, R, C, B_, SIZE(B_, 1), X, SIZE(X, 1), rcond, ferr, berr, work, &
                xwork, info)
    $ASSERT(info == 0,Non-zero return from XGESVX)

    ! Finish

    return

  end function linear_solve_mat_${INFIX}_

  $endsub

  $LINEAR_SOLVE(r,real)
  $LINEAR_SOLVE(c,complex)

  !****

  $define $COMMUTATOR $sub

  $local $INFIX $1
  $local $TYPE $2

  function commutator_${INFIX}_ (A, B) result (C)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: B(:,:)
    $TYPE(WP)             :: C(SIZE(A, 1),SIZE(A, 2))
    
    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))

    $CHECK_BOUNDS(SIZE(B, 1),SIZE(B, 2))
    $CHECK_BOUNDS(SIZE(B, 1),SIZE(A, 1))

    ! Calculate the commutator [A,B]

    C = MATMUL(A, B) - MATMUL(B, A)

    ! Finish

    return

  end function commutator_${INFIX}_

  $endsub

  $COMMUTATOR(r,real)
  $COMMUTATOR(c,complex)

  !****

  $define $MEASURE_BANDWIDTH_FULL $sub
  
  $local $INFIX $1
  $local $TYPE $2
  $local $ZERO $3

  subroutine measure_bandwidth_full_${INFIX}_ (A, n_l, n_u)

    $TYPE(WP), intent(in) :: A(:,:)
    integer, intent(out)  :: n_l
    integer, intent(out)  :: n_u

    integer :: n
    integer :: i
    integer :: j
    integer :: i_min
    integer :: i_max
    integer :: k

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))

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

  end subroutine measure_bandwidth_full_${INFIX}_

  $endsub

  $MEASURE_BANDWIDTH_FULL(r,real,0._WP)
  $MEASURE_BANDWIDTH_FULL(c,complex,(0._WP,0._WP))

  !****

  $define $DIAGONAL $sub

  $local $INFIX $1
  $local $TYPE $2

  pure function diagonal_${INFIX}_ (A) result (D)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP)             :: d(MIN(SIZE(A, 1),SIZE(A, 2)))

    integer :: i

    ! Extract the diagonal elements of matrix A

    do i = 1,SIZE(d)
       d(i) = A(i,i)
    end do

    ! Finish

    return

  end function diagonal_${INFIX}_

  $endsub

  $DIAGONAL(r,real)
  $DIAGONAL(c,complex)

  !****

  $define $DIAGONAL_MATRIX $sub

  $local $INFIX $1
  $local $TYPE $2

  pure function diagonal_matrix_${INFIX}_ (x) result (D)

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

  end function diagonal_matrix_${INFIX}_

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

  $local $INFIX $1
  $local $TYPE $2

  function matrix_exp_${INFIX}_ (A, t) result (exp_At)

    $TYPE(WP), intent(in) :: A(:,:)
    real(WP), intent(in)  :: t
    $TYPE(WP)             :: exp_At(SIZE(A, 1),SIZE(A, 2))

    integer, parameter :: IDEG = 6

    $TYPE(WP)   :: A_(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: lambda(SIZE(A, 1))
    complex(WP) :: V_l(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_r(SIZE(A, 1),SIZE(A, 2))

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))

    ! Calculate the matrix exponential exp(A*t)

    A_ = A

    call eigen_decompose_${INFIX}_(A_, lambda, V_l=V_l, V_r=V_r)

    $if($TYPE eq 'real')
    exp_At = REAL(matrix_exp_eigen_(lambda, V_l, V_r, t))
    $else
    exp_At = matrix_exp_eigen_(lambda, V_l, V_r, t)
    $endif

    ! Finish

    return

  end function matrix_exp_${INFIX}_

  $endsub

  $MATRIX_EXP(r,real)
  $MATRIX_EXP(c,complex)

  !****

  function matrix_exp_eigen_ (lambda, V_l, V_r, t) result (exp_At)

    complex(WP), intent(in) :: lambda(:)
    complex(WP), intent(in) :: V_l(:,:)
    complex(WP), intent(in) :: V_r(:,:)
    real(WP), intent(in)    :: t
    complex(WP)             :: exp_At(SIZE(lambda),SIZE(lambda))

    integer :: i

    $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(V_l, 1))
    $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(V_r, 1))

    $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(lambda))
    $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(lambda))

    ! Calculate the matrix exponential exp(At) from the
    ! eigendecomposition of A

    do i = 1,SIZE(lambda)
       exp_At(i,:) = EXP(lambda(i)*t)*V_l(i,:)
    end do

    exp_At = MATMUL(V_r, exp_At)

    ! Finish

    return

  end function matrix_exp_eigen_

  !****

  $define $OUTER_PRODUCT $sub

  $local $INFIX $1
  $local $TYPE $2

  function outer_product_${INFIX}_ (v, w) result (A)

    $TYPE(WP), intent(in) :: v(:)
    $TYPE(WP), intent(in) :: w(:)
    $TYPE(WP)             :: A(SIZE(v), SIZE(w))

    ! Calculate the outer product v w^T

    A = SPREAD(v, 2, SIZE(w))*SPREAD(w, 1, SIZE(v))

    ! Finish

    return

  end function outer_product_${INFIX}_

  $endsub

  $OUTER_PRODUCT(r,real)
  $OUTER_PRODUCT(c,complex)

  !****

  function outer_product_l_ (v, w) result (A)

    logical, intent(in) :: v(:)
    logical, intent(in) :: w(:)
    logical             :: A(SIZE(v), SIZE(w))

    ! Calculate the outer product v w^T

    A = SPREAD(v, 2, SIZE(w)) .AND. SPREAD(w, 1, SIZE(v))

    ! Finish

    return

  end function outer_product_l_

  !****

  $define $PARTITION $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine partition_${INFIX}_ (A, mask_r, mask_c, A_part)

    $TYPE(WP), intent(in)  :: A(:,:)
    logical, intent(in)    :: mask_r(:)
    logical, intent(in)    :: mask_c(:)
    $TYPE(WP), intent(out) :: A_part(:,:)

    $CHECK_BOUNDS(SIZE(mask_r),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(mask_c),SIZE(A, 2))

    $CHECK_BOUNDS(COUNT(mask_r),SIZE(A_part, 1))
    $CHECK_BOUNDS(COUNT(mask_c),SIZE(A_part, 2))

    ! Partition A according to the row and column masks

    A_part = RESHAPE(PACK(A, outer_product(mask_r, mask_c)), SHAPE(A_part))

    ! Finish

    return

  end subroutine partition_${INFIX}_

  $endsub

  $PARTITION(r,real)
  $PARTITION(c,complex)

end module gyre_linalg
