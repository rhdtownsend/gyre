! Program  : gyre_tar_eigen
! Purpose  : traditional approximation of rotation (TAR) eigenvalue calculation
!
! Copyright 2016 Rich Townsend
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

module gyre_tar_eigen

  ! Uses

  use core_kinds

  use f95_lapack

  use ISO_FORTRAN_ENV

  $if ($IEEE)
  use IEEE_ARITHMETIC
  $endif

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: lambda

  ! Procedures

contains

  function lambda (nu, m, k)

    real(WP), intent(in)     :: nu
    integer, intent(in)      :: m
    integer, intent(in)      :: k
    real(WP)                 :: lambda

    ! Calculate the (m, k) Hough eigenvalue lambda

    lambda = lambda_matrix_(nu, m, k, 'BISECT')

    ! Finish

    return

  end function lambda

  !****

  function lambda_matrix_ (nu, m, k, algo) result (lambda)

    real(WP), intent(in)     :: nu
    integer, intent(in)      :: m
    integer, intent(in)      :: k
    character(*), intent(in) :: algo
    real(WP)                 :: lambda

    logical  :: parity
    integer  :: n
    real(WP) :: lambda_A
    real(WP) :: lambda_A_cmp

    ! Calculate the (m, k) Hough eigenvalue lambda, the using the
    ! matrix approach described in Townsend (1997, PhD Thesis,
    ! University of London) and references therein. The matrix
    ! dimension is determined dynamically

    ! Special-case handling of m=0, k=0

    if (m == 0 .AND. k == 0) then
       lambda = 0._WP
       return
    endif

    ! Initialize n to the smallest possible value

    parity = MOD(k, 2) == 0

    if (k >= 0) then

       if (parity) then
          n = k/2+1
       else
          n = (k+1)/2+1
       endif

       if (m*nu < 0._WP) then
          n = n + Xi_(nu, m, n, parity)
       endif

    else

       if (parity) then
          n = ABS(k)/2
       else
          n = ABS(k-1)/2
       endif

    endif

    ! Double n until a good A-matrix eigenvalue (by convention,
    ! non-zero) is found

    dble_loop : do

       lambda_A = lambda_A_(nu, m, k, n, algo)

       if (lambda_A /= 0._WP) exit dble_loop

       n = 2*n

    end do dble_loop

    ! Continue to double n until convergence is achieved

    conv_loop : do

       lambda_A_cmp = lambda_A_(nu, m, k, 2*n, algo)

       if (lambda_A_cmp /= 0._WP) then

          if (lambda_A == lambda_A_cmp) exit conv_loop

          lambda_A = lambda_A_cmp

       endif

       n = 2*n

    end do conv_loop

    ! Set up lambda

    lambda = 1._WP/lambda_A

    ! Finish

    return

  end function lambda_matrix_

  !****

  function lambda_A_ (nu, m, k, n, algo) result (lambda)

    real(WP), intent(in)     :: nu
    integer, intent(in)      :: m
    integer, intent(in)      :: k
    integer, intent(in)      :: n
    character(*), intent(in) :: algo
    real(WP)                 :: lambda

    ! Calculate the index-k eigenvalue of the A-matrix with fixed
    ! dimension n, using the indicated algorithm

    select case (algo)
    case('LAPACK')
       lambda = lambda_A_lapack_(nu, m, k, n)
    case ('BISECT')
       lambda = lambda_A_bisect_(nu, m, k, n)
    case default
       $ABORT(Unrecognized algo)
    end select

    ! Finish

    return

  end function lambda_A_

  !****

  function lambda_A_lapack_ (nu, m, k, n) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    integer, intent(in)  :: n
    real(WP)             :: lambda

    $if ($IEEE)
    logical  :: h_over
    logical  :: h_zero
    logical  :: h_inval
    $endif
    logical  :: parity
    real(WP) :: A_1(n)
    real(WP) :: A_2n(n)
    real(WP) :: A_2d(n)
    real(WP) :: A_3n(n)
    real(WP) :: A_3d(n)
    real(WP) :: A_4n(n)
    real(WP) :: A_4d(n)
    real(WP) :: A_5(n)
    real(WP) :: A_D(n)
    real(WP) :: A_E(n)
    integer  :: i
    real(WP) :: W(n)
    integer  :: j
    integer  :: j_inf
    
    ! Calculate the index-k eigenvalue of the A-matrix with fixed
    ! dimension n, using LAPACK routines at working-precision

    $if ($IEEE)
    call save_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Determine the (negative-to-positive) index of the eigenvalue to
    ! calculate

    parity = MOD(k, 2) == 0

    if (k >= 0) then
       if (parity) then
          i = n - k/2
       else
          i = n - (k-1)/2
       endif
    else
       if (parity) then
          i = ABS(k)/2
       else
          i = ABS(k-1)/2
       endif
    endif

    i = MODULO(i-Xi_(nu, m, n, parity)-1, n) + 1

    ! Assemble the A-matrix

    call assemble_A_WP_(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, &
                        A_D, A_E)
    
    A_E(n) = 0._WP

    ! Calculate the eigenvalue, taking special care to handle cases
    ! where elements of A are infinite

    if (A_5(1) == 0._WP) then

       ! A(1,1) infinite, causing a +inf eigenvalue

       if (n > 1 .AND. i < n) then

          call LA_STEVR(A_D(2:), A_E(2:), W(:n-1), IL=i, IU=i)
          lambda = W(1)

       else

          $if ($IEEE)
          lambda = IEEE_VALUE(0._WP, IEEE_POSITIVE_INF)
          $else
          lambda = HUGE(0._WP)
          $endif

       endif

    elseif (A_2d(1) == 0._WP .AND. A_2n(1) /= 0._WP) then

       ! A(1,1) infinite, causing a -inf eigenvalue

       if (n > 1 .AND. i > 1) then

          call LA_STEVR(A_D(2:), A_E(2:), W(2:), IL=i-1, IU=i-1)
          lambda = W(2)

       else

          $if ($IEEE)
          lambda = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)
          $else
          lambda = -HUGE(0._WP)
          $endif

       endif

    elseif (A_3d(n) == 0._WP .AND. A_3n(n) /= 0._WP) then

       ! A(n,n) infinite, causing a -inf eigenvalue

       if (n > 1 .AND. i > 1) then

          call LA_STEVR(A_D(:n-1), A_E(:n-1), W(2:), IL=i-1, IU=i-1)
          lambda = W(2)

       else

          $if ($IEEE)
          lambda = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)
          $else
          lambda = -HUGE(0._WP)
          $endif

       endif

    else
       
       ! Scan for a 2x2 block of infinite components (arising from a
       ! pure Rossby mode being part of the spectrum)

       j_inf = 0

       scan_loop : do j = 1, n-1
          if (A_3d(j) == 0._WP .AND. A_3n(j) /= 0._WP) then
             j_inf = j
             exit scan_loop
          endif
       end do scan_loop

       if (j_inf /= 0) then

          print *,nu,m,k,j_inf

          $ABORT(Pure Rossby mode encountered)

       else

          call LA_STEVR(A_D, A_E, W, IL=i, IU=i)
          lambda = W(1)

       endif

    endif

    $if ($IEEE)
    call restore_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Finish

    return

  end function lambda_A_lapack_

  !****

  function lambda_A_bisect_ (nu, m, k, n) result (lambda)

    real(WP), intent(in)  :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    integer, intent(in)   :: n
    real(WP)              :: lambda

    real(WP), parameter :: NU_TRANS = 1E4_WP

    $if ($IEEE)
    logical  :: h_over
    logical  :: h_zero
    logical  :: h_inval
    $endif
    real(WP) :: lambda_est
    logical  :: parity
    real(QP) :: A_1(n)
    real(QP) :: A_2n(n)
    real(QP) :: A_2d(n)
    real(QP) :: A_3n(n)
    real(QP) :: A_3d(n)
    real(QP) :: A_4n(n)
    real(QP) :: A_4d(n)
    real(QP) :: A_5(n)
    real(QP) :: A_D(n)
    real(QP) :: A_E(n)
    
    ! Calculate the index-k eigenvalue of the A-matrix with fixed
    ! dimension n, using bisection routines at quad-precision

    $if ($IEEE)
    call save_ieee_state_(h_over, h_zero, h_inval)
    $endif
 
    ! Get an initial estimate for the eigenvalue; use the LAPACK
    ! routine for |nu| < NU_TRANS, and the asymptotic value plus an
    ! extrapolated error correction otherwise

    if (ABS(nu) < NU_TRANS) then
       lambda_est = lambda_A_lapack_(nu, m, k, n)
    else
       lambda_est = 1._WP/(lambda_asymp_(nu, m, k) + &
                           (1._WP/lambda_A_lapack_(SIGN(NU_TRANS, nu), m, k, n) - lambda_asymp_(SIGN(NU_TRANS, nu), m, k))/nu**2)
    endif
    
    ! Check to see if lambda_est is +/- HUGE, indicating a Rossby mode

    if (lambda_est == HUGE(0._WP) .OR. lambda_est == -HUGE(0._WP)) then
       lambda = lambda_est
       return
    end if

    ! Assemble the A-matrix

    parity = MOD(k, 2) == 0

    call assemble_A_QP_(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, &
                        A_D, A_E)

    ! Calculate the eigenvalue, taking special care to handle cases
    ! where elements of A are infinite

    if (A_5(1) == 0._QP) then

       ! A(1,1) infinite, causing a +inf eigenvalue

       lambda = bisect_lambda_(A_D(2:), A_E(2:), lambda_est)

    elseif (A_2d(1) == 0._QP .AND. A_2n(1) /= 0._QP) then

       ! A(1,1) infinite, causing a -inf eigenvalue

       lambda = bisect_lambda_(A_D(2:), A_E(2:), lambda_est)

    elseif (A_3d(n) == 0._WP .AND. A_3n(n) /= 0._WP) then

       ! A(n,n) infinite, causing a -inf eigenvalue

       lambda = bisect_lambda_(A_D(:n-1), A_E(:n-1), lambda_est)

    else
       
       lambda = bisect_lambda_(A_D, A_E, lambda_est)

    endif

    ! Finish

    return

  contains

    function bisect_lambda_ (A_D, A_E, lambda_est) result (lambda)

      real(QP), intent(in) :: A_D(:)
      real(QP), intent(in) :: A_E(:)
      real(WP), intent(in) :: lambda_est
      real(WP)             :: lambda

      real(QP) :: lambda_a
      real(QP) :: lambda_b
      real(QP) :: lambda_m
      real(QP) :: dlambda
      integer  :: P_a
      integer  :: P_b
      integer  :: P_m
      
      $CHECK_BOUNDS(SIZE(A_E),SIZE(A_D))

      ! Set up a bracket around the estimated eigenvalue lambda_est

      if (lambda_est > 0._WP) then
         lambda_a = lambda_est*(1._WP - SQRT(EPSILON(0._WP)))
         lambda_b = lambda_est*(1._WP + SQRT(EPSILON(0._WP)))
      else
         lambda_a = lambda_est*(1._WP + SQRT(EPSILON(0._WP)))
         lambda_b = lambda_est*(1._WP - SQRT(EPSILON(0._WP)))
      endif

      expand_loop : do

         P_a = bisect_P_(A_D, A_E, lambda_a)
         P_b = bisect_P_(A_D, A_E, lambda_b)

         if (P_b == P_a+1) then

            ! Single eigenvalue found

            exit expand_loop

         elseif (P_b > P_a+1) then

            ! Multiple eigenvalues found --- return a bogus value
            ! (which has the ultimate effect that we'll try again at
            ! larger n)

            lambda = 0._WP
            return

         endif

         ! Expand the bracket

         lambda_m = 0.5_QP*(lambda_a + lambda_b)
         dlambda = lambda_b - lambda_a

         lambda_a = lambda_a - 0.5_WP*dlambda
         lambda_b = lambda_b + 0.5_WP*dlambda

      end do expand_loop

      ! Now apply the bisection algorithm

      bisect_loop : do

         ! Check for convergence

         if (REAL(lambda_a, WP) == REAL(lambda_b, WP)) then
            lambda = REAL(lambda_a, WP)
            exit bisect_loop
         endif

         ! Set up the midpoint

         lambda_m = 0.5_QP*(lambda_a + lambda_b)
         
         ! Calculate P at the midpoint
       
         P_m = bisect_P_(A_D, A_E, lambda_m)

         ! Update the bracket
         
         if (P_m == P_a) then
            lambda_a = lambda_m
         elseif (P_m == P_b) then
            lambda_b = lambda_m
         else
            print *,P_a, P_b, P_m
            $ABORT(P_m not equal to P_a or P_b)
         endif
         
      end do bisect_loop

      ! Finish

      return

    end function bisect_lambda_

    function bisect_P_ (A_D, A_E, lambda) result (P)

      real(QP), intent(in) :: A_D(:)
      real(QP), intent(in) :: A_E(:)
      real(QP), intent(in) :: lambda
      integer              :: P

      real(QP), parameter :: RELFEH = EPSILON(0._QP)

      real(QP) :: q(SIZE(A_D))
      integer  :: i

      $CHECK_BOUNDS(SIZE(A_E),SIZE(A_D))

      $ASSERT(SIZE(A_D) > 0,Empty matrix)

      ! Calculate the number P of eigenvalues of A which are less than
      ! lambda, using the algorithm in Barth, Martin & Wilkinson
      ! (1967)

      q(1) = A_D(1) - lambda

      do i = 2, SIZE(A_D)
         if (q(i-1) /= 0._QP) then
            q(i) = (A_D(i) - lambda) - A_E(i-1)**2/q(i-1)
         else
            q(i) = (A_D(i) - lambda) - ABS(A_E(i-1))/RELFEH
         endif
      end do

      P = COUNT(q < 0._QP)

      ! Finish

      return

    end function bisect_P_

  end function lambda_A_bisect_

  !****

  function lambda_asymp_ (nu, m, k) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda

    integer  :: s

    ! Calculate the (m, k) Hough eigenvalue lambda, using asymptotic expansions

    if (m*nu >= 0._WP) then

       ! Prograde modes

       if (k >= 0) then

          s = k - 1

          if (s >= 0) then
             lambda = lambda_asymp_grav_(nu, m, s)
          else
             lambda = lambda_asymp_kelv_(nu, m)
          endif

       else

          $ABORT(Invalid k for prograde modes)

       endif

    else

       ! Retrograde modes

       if (k >= -1) then

          s = k + 1

          lambda = lambda_asymp_grav_(nu, m, s)

       else

          s = - k - 1

          lambda = lambda_asymp_ross_(nu, m, s)

       endif

    end if

  contains

    function lambda_asymp_grav_ (nu, m, s) result (lambda)

      real(WP), intent(in) :: nu
      integer, intent(in)  :: m
      integer, intent(in)  :: s
      real(WP)             :: lambda

      real(QP) :: w
      real(QP) :: I2
      real(QP) :: I4
      real(QP) :: alpha_0
      real(QP) :: alpha_1
      real(QP) :: alpha_2
      real(QP) :: alpha_3
      real(QP) :: alpha

      ! Calculate the (m, s) gravity-wave Hough eigenvalue lambda,
      ! using the third-order perturbation expansion by Townsend (in
      ! prep)

      w = 1._QP/REAL(nu, QP)
    
      I2 = 0.5_QP*(2*s+1)
      I4 = 0.75_QP*(2*s*(s+1) + 1)

      alpha_0 = 2*s + 1
      alpha_1 = m/alpha_0
      alpha_2 = m**2/alpha_0 - m*alpha_1/alpha_0**2 + 2/alpha_0*(I4 - alpha_0*I2) + I2 - I4/alpha_0
      alpha_3 = m*(alpha_1**2 - alpha_0*alpha_2)/alpha_0**3 - m**2*alpha_1/alpha_0**2 - &
                2*alpha_1/alpha_0**2*(I4 - alpha_0*I2) -  m/alpha_0**2*I2 + alpha_1/alpha_0**2*I4

      alpha = alpha_0 + alpha_1*w + alpha_2*w**2 + alpha_3*w**3

      lambda = REAL(alpha**2/w**2, WP)

      ! Finish

      return

    end function lambda_asymp_grav_

    function lambda_asymp_kelv_ (nu, m) result (lambda)

      real(WP), intent(in) :: nu
      integer, intent(in)  :: m
      real(WP)             :: lambda

      real(QP) :: w
      real(QP) :: alpha_0
      real(QP) :: alpha_1
      real(QP) :: alpha

      ! Calculate the (m) Kelvin-wave Hough eigenvalue lambda, using
      ! the first-order perturbation expansion by Townsend (in prep)

      alpha_0 = ABS(m)
      alpha_1 = -0.25_QP

      alpha = alpha_0 + alpha_1*w

      lambda = REAL(alpha**2, WP)
      
      ! Finish

      return

    end function lambda_asymp_kelv_

    function lambda_asymp_ross_ (nu, m, s) result (lambda)

      real(WP), intent(in) :: nu
      integer, intent(in)  :: m
      integer, intent(in)  :: s
      real(WP)             :: lambda

      real(QP) :: w
      real(QP) :: alpha_0
      real(QP) :: alpha_1
      real(QP) :: alpha

      ! Calculate the (m, s) Rossby-wave Hough eigenvalue lambda,
      ! using the first-order perturbation expansion by Townsend (in
      ! prep)

      w = 1._QP/REAL(nu, QP)

      alpha_0 = REAL(m, QP)/(2*s+1)
      alpha_1 = REAL(-1 - 4*(1+m**2)*s*(s+1), QP)/(2*s+1)**3

      alpha = alpha_0 + alpha_1*w

      lambda = REAL(alpha**2, WP)

      ! Finish

      return

    end function lambda_asymp_ross_

  end function lambda_asymp_

  !****

  $define $ASSEMBLE_A $sub

  $local $KP $1

  subroutine assemble_A_${KP}_ (nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, A_D, A_E)

    real(WP), intent(in)   :: nu
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    logical, intent(in)    :: parity
    real($KP), intent(out) :: A_1(:)
    real($KP), intent(out) :: A_2n(:)
    real($KP), intent(out) :: A_2d(:)
    real($KP), intent(out) :: A_3n(:)
    real($KP), intent(out) :: A_3d(:)
    real($KP), intent(out) :: A_4n(:)
    real($KP), intent(out) :: A_4d(:)
    real($KP), intent(out) :: A_5(:)
    real($KP), intent(out) :: A_D(:)
    real($KP), intent(out) :: A_E(:)

    real($KP) :: nu_
    integer   :: j
    integer   :: l_j
    real($KP) :: A_2
    real($KP) :: A_3
    real($KP) :: A_4

    $CHECK_BOUNDS(SIZE(A_1),n)
    $CHECK_BOUNDS(SIZE(A_2n),n)
    $CHECK_BOUNDS(SIZE(A_2d),n)
    $CHECK_BOUNDS(SIZE(A_3n),n)
    $CHECK_BOUNDS(SIZE(A_3d),n)
    $CHECK_BOUNDS(SIZE(A_4n),n)
    $CHECK_BOUNDS(SIZE(A_4d),n)
    $CHECK_BOUNDS(SIZE(A_5),n)

    $CHECK_BOUNDS(SIZE(A_D),n)
    $CHECK_BOUNDS(SIZE(A_E),n)

    ! Evaluate constituent parts of the A matrix. These are defined so
    ! that:
    !
    ! A_D(j) = A(j,j) = A_D(j) = [A_1(j) + A_2n(j)/A_2d(j) + A_3n(j)/A_3d(j)] / A_5(j)
    ! A_E(j) = A(j+1,j) = A(j, j+1) = A_4n(j)/A_4d(j)
    !
    ! A_E(n) is unused

    nu_ = REAL(nu, ${KP})

    if (parity) then
       l_j = ABS(m)
    else
       l_j = ABS(m) + 1
    endif
        
    do j = 1, n

       if (l_j == ABS(m)) then
          A_1(j) = 1._${KP} + SIGN(1, m)*nu_/REAL(l_j+1, ${KP})
          A_2n(j) = 0._${KP}
       else
          A_1(j) = 1._${KP} + m*nu_/(REAL(l_j, ${KP})*REAL(l_j+1, ${KP}))
          A_2n(j) = -nu_**2*REAL(l_j-1, ${KP})**2*REAL(l_j+1, ${KP})*J_lm_${KP}_(l_j, m)**2/REAL(l_j, ${KP})
       endif

       A_2d(j) = REAL(l_j-1, ${KP})*REAL(l_j, ${KP}) + m*nu_
       
       A_3n(j) = -nu_**2*REAL(l_j, ${KP})*REAL(l_j+2, ${KP})**2*J_lm_${KP}_(l_j+1, m)**2/REAL(l_j+1, ${KP})
       A_3d(j) = REAL(l_j+1, ${KP})*REAL(l_j+2, ${KP}) + m*nu_

       A_4n(j) = -nu_**2*J_lm_${KP}_(l_j+1, m)*J_lm_${KP}_(l_j+2, m)
       A_4d(j) = A_3d(j)

       A_5(j) = REAL(l_j, ${KP})*REAL(l_j+1, ${KP})

       l_j = l_j + 2

    end do

    ! Now assemble the matrix

    do j = 1, n-1

       if (A_2n(j) /= 0._${KP}) then
          A_2 = A_2n(j)/A_2d(j)
       else
          A_2 = 0._${KP}
       endif

       if (A_3n(j) /= 0._${KP}) then
          A_3 = A_3n(j)/A_3d(j)
       else
          A_3 = 0._${KP}
       endif

       if (A_4n(j) /= 0._${KP}) then
          A_4 = A_4n(j)/A_4d(j)
       else
          A_4 = 0._${KP}
       endif

       if (A_5(j) /= 0._${KP}) then
          A_D(j) = (A_1(j) + A_2 + A_3)/A_5(j)
       endif
       A_E(j) = A_4

    end do

    if (A_2n(n) /= 0._${KP}) then
       A_2 = A_2n(n)/A_2d(n)
    else
       A_2 = 0._${KP}
    endif

    if (A_3n(n) /= 0._${KP} .AND. parity) then
       A_3 = A_3n(n)/A_3d(n)
    else
       A_3 = 0._${KP}
    endif

    if (A_5(n) /= 0._${KP}) then
       A_D(n) = (A_1(n) + A_2 + A_3)/A_5(n)
    endif

    ! Finish

    return

  end subroutine assemble_A_${KP}_

  $endsub

  $ASSEMBLE_A(WP)
  $ASSEMBLE_A(QP)

  !****

  $define $J_LM $sub

  $local $KP $1

  function J_lm_${KP}_ (l, m) result (J_lm_)

    integer, intent(in) :: l
    integer, intent(in) :: m
    real($KP)           :: J_lm_

    ! Calculate the J_lm function, being the integral of Y_lm Y_{l+1}m
    ! cos(theta) over all solid angles

    if (ABS(m) >= l) then

       J_lm_ = 0._${KP}

    else

       J_lm_ = SQRT(REAL(l-m, ${KP})*REAL(l+m, ${KP})/(REAL(2*l+1, ${KP})*REAL(2*l-1, ${KP})))

    endif

    ! Finish

    return
    
  end function J_lm_${KP}_

  $endsub

  $J_LM(WP)
  $J_LM(QP)

  !****

  function Xi_ (nu, m, n, parity)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    logical, intent(in)  :: parity
    integer              :: Xi_

    integer :: l_j
    integer :: j

    ! Calculate the permutation constant Xi (see Townsend 1997)

    Xi_ = 0

    if (-m*nu > 0._WP) then

       if (parity) then
          l_j = ABS(m)
       else
          l_j = ABS(m) + 1
       endif
       
       j_loop : do j = 1,n
          if (-m*nu > REAL(l_j+1, WP)*REAL(l_j+2, WP)) Xi_ = Xi_ + 1
          l_j = l_j + 2
       enddo j_loop
       
       if (.NOT. parity .AND. -m*nu > REAL(ABS(m), WP)*REAL(ABS(m)+1, WP)) Xi_ = Xi_ + 1
       
    endif

    ! Finish

    return

  end function Xi_

  !****

  $if ($IEEE)

  subroutine save_ieee_state_ (h_over, h_zero, h_inval)

    logical, intent(out) :: h_over
    logical, intent(out) :: h_zero
    logical, intent(out) :: h_inval

    ! Save the IEEE state, and then turn off exceptions

    call IEEE_GET_HALTING_MODE(IEEE_OVERFLOW, h_over)
    call IEEE_GET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, h_zero)
    call IEEE_GET_HALTING_MODE(IEEE_INVALID, h_inval)

    call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE.)
    call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE.)
    call IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)

    ! Finish

    return

  end subroutine save_ieee_state_

  !****

  subroutine restore_ieee_state_ (h_over, h_zero, h_inval)

    logical, intent(in) :: h_over
    logical, intent(in) :: h_zero
    logical, intent(in) :: h_inval

    ! Restore the IEEE state

    call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, h_over)
    call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, h_zero)
    call IEEE_SET_HALTING_MODE(IEEE_INVALID, h_inval)

    ! Finish

    return

  end subroutine restore_ieee_state_

  $endif

end module gyre_tar_eigen
