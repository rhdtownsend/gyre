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
     module procedure eigen_decompose_r
     module procedure eigen_decompose_c
  end interface eigen_decompose

  interface linear_solve
     module procedure linear_solve_vec_r
     module procedure linear_solve_vec_c
     module procedure linear_solve_mat_r
     module procedure linear_solve_mat_c
  end interface linear_solve

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

  public :: eigen_decompose
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
    real(WP)    :: work(4*SIZE(A,1))
    integer     :: info
    complex(WP) :: V_l_c(SIZE(A, 1),SIZE(A, 2))
    complex(WP) :: V_r_c(SIZE(A, 1),SIZE(A, 2))
    integer     :: i
    integer     :: ind(SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(lambda),SIZE(A, 1))

    if(PRESENT(V_l)) then
       $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(A, 2))
    endif
    
    if(PRESENT(V_r)) then
       $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(A, 2))
    endif

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    A_ = A

    if(PRESENT(V_l) .OR. PRESENT(V_r)) then

       call XGEEV('V', 'V', SIZE(A, 1), A_, SIZE(A, 1), lambda_re, lambda_im, &
                  V_l_r, SIZE(V_l_r, 1), V_r_r, SIZE(V_r_r, 1), &
                  work, SIZE(work), info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

       call reconstruct_eigenvector(lambda_re, lambda_im, V_l_r, V_l_c)
       call reconstruct_eigenvector(lambda_re, lambda_im, V_r_r, V_r_c)

       V_l_c = CONJG(TRANSPOSE(V_l_c))

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       do i = 1,SIZE(A,1)
          V_l_c(i,:) = V_l_c(i,:)/SUM(V_l_c(i,:)*V_r_c(:,i))
       enddo

       if(PRESENT(V_l)) V_l = V_l_c
       if(PRESENT(V_r)) V_r = V_r_c

    else

       call XGEEV('N', 'N', SIZE(A, 1), A_, SIZE(A, 1), lambda_re, lambda_im, &
                  V_l_r, SIZE(V_l_r, 1), V_r_r, SIZE(V_r_r, 1), &
                  work, SIZE(work), info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

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
    complex(WP) :: work(2*SIZE(A, 1))
    real(WP)    :: rwork(2*SIZE(A, 1))
    integer     :: info
    integer     :: i
    integer     :: ind(SIZE(A, 1))

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(lambda),SIZE(A, 1))

    if(PRESENT(V_l)) then
       $CHECK_BOUNDS(SIZE(V_l, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_l, 2),SIZE(A, 2))
    endif
    
    if(PRESENT(V_r)) then
       $CHECK_BOUNDS(SIZE(V_r, 1),SIZE(A, 1))
       $CHECK_BOUNDS(SIZE(V_r, 2),SIZE(A, 2))
    endif

    if(PRESENT(sort)) then
       sort_ = sort
    else
       sort_ = .FALSE.
    end if

    ! Perform the eigendecomposition of A

    A_ = A

    if(PRESENT(V_l) .OR. PRESENT(V_r)) then

       call XGEEV('V', 'V', SIZE(A, 1), A_, SIZE(A, 1), lambda, &
                  V_l_c, SIZE(V_l_c, 1), V_r_c, SIZE(V_r_c, 1), &
                  work, SIZE(work, 1), rwork, info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

       V_l_c = CONJG(TRANSPOSE(V_l_c))

       ! Renormalize the left eigenvectors so they are orthonormal to
       ! the right eigenvectors

       do i = 1,SIZE(A, 1)
          V_l_c(i,:) = V_l_c(i,:)/SUM(V_l_c(i,:)*V_r_c(:,i))
       enddo

       if(PRESENT(V_r)) V_r = V_r_c
       if(PRESENT(V_l)) V_l = V_l_c
       
    else

       call XGEEV('N', 'N', SIZE(A, 1), A_, SIZE(A, 1), lambda, &
                  V_l_c, SIZE(V_l_c, 1), V_r_c, SIZE(V_r_c, 1), &
                  work, SIZE(work, 1), rwork, info)
       $ASSERT(info == 0,Non-zero return from XGEEV)

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

    $TYPE(WP) :: Mx(SIZE(x), 1)
    $TYPE(WP) :: Mb(SIZE(b), 1)

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(b),SIZE(A, 1))

    ! Solve the linear system A x = b

    Mx = linear_solve(A, Mb)

    x = Mx(:,1)

    ! Finish

    return

  end function linear_solve_vec_$SUFFIX

  function linear_solve_mat_$SUFFIX (A, B) result (X)

    $TYPE(WP), intent(in) :: A(:,:)
    $TYPE(WP), intent(in) :: B(:,:)
    $TYPE(WP)             :: X(SIZE(A, 2), SIZE(B, 2))

    $TYPE(WP) :: A_(SIZE(A, 1),SIZE(A, 2))
    integer   :: ipiv(SIZE(A, 1))
    integer   :: info

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))
    $CHECK_BOUNDS(SIZE(B, 1),SIZE(A, 1))

    ! Solve the linear system A X = B

    A_ = A
    X = B

    call XGESV(SIZE(A, 1), SIZE(B, 2), A_, SIZE(A, 1), ipiv, X, SIZE(B, 1), info)
    $ASSERT(info == 0,Non-zero return from XGESV)

    ! Finish

    return

  end function linear_solve_mat_$SUFFIX

  $endsub

  $LINEAR_SOLVE(r,real)
  $LINEAR_SOLVE(c,complex)

!****

  $define $COMMUTATOR $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function commutator_$SUFFIX (A, B) result (C)

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

    $CHECK_BOUNDS(SIZE(A, 1),SIZE(A, 2))

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

    $CHECK_BOUNDS(SIZE(mask_r),SIZE(A, 1))
    $CHECK_BOUNDS(SIZE(mask_c),SIZE(A, 2))

    $CHECK_BOUNDS(COUNT(mask_r),SIZE(A_part, 1))
    $CHECK_BOUNDS(COUNT(mask_c),SIZE(A_part, 2))

    ! Partition A according to the row and column masks

    A_part = RESHAPE(PACK(A, outer_product(mask_r, mask_c)), SHAPE(A_part))

    ! Finish

    return

  end subroutine partition_$SUFFIX

  $endsub

  $PARTITION(r,real)
  $PARTITION(c,complex)

end module gyre_linalg
