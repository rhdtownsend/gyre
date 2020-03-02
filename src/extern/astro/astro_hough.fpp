! Module  : astro_hough
! Purpose : calculation of Hough eigenvalues & eigenfunctions, incl. hough_t type
! Notes   : This code adopts the sign convention m*nu > 0 for prograde modes

$include 'core.inc'

module astro_hough

  ! Uses

  use core_kinds
  use core_linalg
  use core_order

  use f95_lapack

  use ISO_FORTRAN_ENV

  $if ($IEEE)
  use IEEE_ARITHMETIC
  $endif
  
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: hough_t
     private
     real(WP), allocatable :: b(:)
     real(WP), allocatable :: c_r(:)
     real(WP), allocatable :: c_t(:)
     real(WP), allocatable :: c_p(:)
     real(WP), public      :: lambda
     real(WP), public      :: nu
     integer, public       :: m
     integer, public       :: k
     integer               :: n
   contains
     private
     procedure, public :: Theta => Theta_
     procedure, public :: Theta_xi => Theta_xi_
  end type hough_t

  ! Interfaces

  interface hough_t
     module procedure hough_t_n_
     module procedure hough_t_tol_
  end interface hough_t

  interface lambda
     module procedure lambda_n_
     module procedure lambda_tol_
  end interface lambda

  ! Access specifiers

  private

  public :: hough_t
  public :: lambda
  public :: lambda_asymp

  ! Procedures

contains

  function hough_t_n_ (nu, m, k, n) result (ho)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    integer, intent(in)  :: n
    type(hough_t)        :: ho

    ! Construct the hough_t for fixed matrix dimension

    call eigen_W_n_(nu, m, k, n, ho%lambda, ho%b)

    allocate(ho%c_r(n))
    allocate(ho%c_t(n))
    allocate(ho%c_p(n))

    call eval_coeffs_(nu, m, k, ho%lambda, ho%b, ho%c_r, ho%c_t, ho%c_p)

    ho%nu = nu

    ho%m = m
    ho%k = k
    ho%n = n

    ! Finish

    return

  end function hough_t_n_

!****

  function hough_t_tol_ (nu, m, k, tol) result (ho)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: tol
    type(hough_t)        :: ho

    integer :: n

    ! Construct the hough_t to the specified tolerance

    call eigen_W_tol_(nu, m, k, tol, ho%lambda, ho%b)

    n = SIZE(ho%b)

    allocate(ho%c_r(n))
    allocate(ho%c_t(n))
    allocate(ho%c_p(n))

    call eval_coeffs_(nu, m, k, ho%lambda, ho%b, ho%c_r, ho%c_t, ho%c_p)

    ho%nu = nu

    ho%m = m
    ho%k = k
    ho%n = n

    ! Finish

    return

  end function hough_t_tol_

!****

  function Theta_ (this, mu) result (Theta)

    class(hough_t), intent(in) :: this
    real(WP), intent(in)       :: mu
    real(WP)                   :: Theta(3)

    logical  :: parity
    real(WP) :: P_lm(2*this%n)

    ! Evaluate the Hough Theta functions (see Townsend 2003) at
    ! mu=cos theta

    parity = MOD(this%k, 2) == 0

    P_lm = P_lm_series(2*this%n, this%m, mu)

    if (parity) then
       associate (P_l_j => P_lm(1::2), P_n_j => P_lm(2::2))
         Theta(1) = DOT_PRODUCT(this%c_r, P_l_j)
         Theta(2) = DOT_PRODUCT(this%c_t, P_n_j)
         Theta(3) = DOT_PRODUCT(this%c_p, P_l_j)
       end associate
    else
       associate (P_l_j => P_lm(2::2), P_n_j => P_lm(1::2))
         Theta(1) = DOT_PRODUCT(this%c_r, P_l_j)
         Theta(2) = DOT_PRODUCT(this%c_t, P_n_j)
         Theta(3) = DOT_PRODUCT(this%c_p, P_l_j)
       end associate
    endif

    ! Finish

    return

  end function Theta_

!****

  function Theta_xi_ (this, mu) result (Theta_xi)

    class(hough_t), intent(in) :: this
    real(WP), intent(in)       :: mu
    real(WP)                   :: Theta_xi(3)

    logical  :: parity
    real(WP) :: stheta
    real(WP) :: P_lm(2*this%n)
    real(WP) :: dP_lm(2*this%n)

    ! Evaluate the Hough theta-function parts of the displacement
    ! vector (see Townsend 2003) at mu=cos theta

    parity = MOD(this%k, 2) == 0

    stheta = SQRT(1._WP - mu**2)

    if (stheta > EPSILON(0._WP)) then

       P_lm = P_lm_series(2*this%n, this%m, mu)
    
       if (parity) then
          associate (P_l_j => P_lm(1::2), P_n_j => P_lm(2::2))
            Theta_xi(1) = DOT_PRODUCT(this%c_r, P_l_j)
            Theta_xi(2) = DOT_PRODUCT(this%c_t, P_n_j)/stheta
            Theta_xi(3) = DOT_PRODUCT(this%c_p, P_l_j)/stheta
          end associate
       else
          associate (P_l_j => P_lm(2::2), P_n_j => P_lm(1::2))
            Theta_xi(1) = DOT_PRODUCT(this%c_r, P_l_j)
            Theta_xi(2) = DOT_PRODUCT(this%c_t, P_n_j)/stheta
            Theta_xi(3) = DOT_PRODUCT(this%c_p, P_l_j)/stheta
          end associate
       endif

    else

       ! Use L'Hopital's rule to deal with the poles

       P_lm = P_lm_series(2*this%n, this%m, mu)
       dP_lm = dP_lm_series(2*this%n, this%m, mu)
       
       if (parity) then
          associate (P_l_j => P_lm(1::2), P_n_j => P_lm(2::2), dP_l_j => dP_lm(1::2), dP_n_j => dP_lm(2::2))
            Theta_xi(1) = DOT_PRODUCT(this%c_r, P_l_j)
            Theta_xi(2) = DOT_PRODUCT(this%c_t, dP_n_j)/(-mu)
            Theta_xi(3) = DOT_PRODUCT(this%c_p, dP_l_j)/(-mu)
          end associate
       else
          associate (P_l_j => P_lm(2::2), P_n_j => P_lm(1::2), dP_l_j => dP_lm(2::2), dP_n_j => dP_lm(1::2))
            Theta_xi(1) = DOT_PRODUCT(this%c_r, P_l_j)
            Theta_xi(2) = DOT_PRODUCT(this%c_t, dP_n_j)/(-mu)
            Theta_xi(3) = DOT_PRODUCT(this%c_p, dP_l_j)/(-mu)
          end associate
       endif

    end if

    ! Finish

    return

  end function Theta_xi_

!****

  subroutine eval_coeffs_ (nu, m, k, lambda, b, c_r, c_t, c_p)

    real(WP), intent(in)  :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    real(WP), intent(in)  :: lambda
    real(WP), intent(in)  :: b(:)
    real(WP), intent(out) :: c_r(:)
    real(WP), intent(out) :: c_t(:)
    real(WP), intent(out) :: c_p(:)

    integer  :: n
    logical  :: parity
    real(WP) :: N_lm(2*SIZE(b))
    real(WP) :: q(SIZE(b))
    real(WP) :: r(SIZE(b))
    integer  :: l_j
    integer  :: n_j
    integer  :: j
    real(WP) :: P_lm(2*SIZE(b))
    real(WP) :: dP_lm(2*SIZE(b))

    $CHECK_BOUNDS(SIZE(c_r),SIZE(b))
    $CHECK_BOUNDS(SIZE(c_t),SIZE(b))
    $CHECK_BOUNDS(SIZE(c_p),SIZE(b))

    ! Calculate the assoc. Legendre polynomial expansion coefficients
    ! for the Hough functions

    n = SIZE(b)

    parity = MOD(k, 2) == 0

    ! Set up the radial coefficients c_r

    N_lm = N_lm_series(2*n, m)

    c_r = N_lm(1::2)*b

    ! Set up the polar (c_t) and azimuthal (c_p) coefficients

    call eval_eigvecs_(nu, m, k, lambda, b, q, r)

    if (parity) then

       l_j = ABS(m)
       n_j = l_j + 1

       associate (N_l_j => N_lm(1::2), N_n_j => N_lm(2::2))

         if (n == 1) then

            c_t(1) = (n_j-1)*(n_j-m)*N_l_j(1)*q(1)/(2*n_j-1) + m*N_n_j(1)*r(1)
            c_p(1) = (l_j+2)*(l_j+m+1)*N_n_j(1)*r(1)/(2*l_j+3) - m*N_l_j(1)*q(1)

         else

            c_t(1) = (n_j-1)*(n_j-m)*N_l_j(1)*q(1)/(2*n_j-1) - (n_j+2)*(n_j+m+1)*N_l_j(2)*q(2)/(2*n_j+3) + &
                     m*N_n_j(1)*r(1)
            c_p(1) = (l_j+2)*(l_j+m+1)*N_n_j(1)*r(1)/(2*l_j+3) - m*N_l_j(1)*q(1)
          
            l_j = l_j + 2
            n_j = n_j + 2

            do j = 2,n-1
               c_t(j) = (n_j-1)*(n_j-m)*N_l_j(j)*q(j)/(2*n_j-1) - (n_j+2)*(n_j+m+1)*N_l_j(j+1)*q(j+1)/(2*n_j+3) + &
                        m*N_n_j(j)*r(j)
               c_p(j) = (l_j+2)*(l_j+m+1)*N_n_j(j)*r(j)/(2*l_j+3) - (l_j-1)*(l_j-m)*N_n_j(j-1)*r(j-1)/(2*l_j-1) - &
                        m*N_l_j(j)*q(j)
               l_j = l_j + 2
               n_j = n_j + 2
            end do

            c_t(n) = (n_j-1)*(n_j-m)*N_l_j(n)*q(n)/(2*n_j-1) - m*N_n_j(n)*r(n)
            c_p(n) = (l_j+2)*(l_j+m+1)*N_n_j(n)*r(n)/(2*l_j+3) - (l_j-1)*(l_j-m)*N_n_j(n-1)*r(n-1)/(2*l_j-1) - &
                     m*N_l_j(n)*q(n)

         endif

       end associate

    else

       n_j = ABS(m)
       l_j = n_j + 1

       associate (N_l_j => N_lm(2::2), N_n_j => N_lm(1::2))

         if (n == 1) then

            c_t(1) = -(n_j+2)*(n_j+m+1)*N_l_j(1)*q(1)/(2*n_j+3) + m*N_n_j(1)*r(1)
            c_p(1) = -(l_j-1)*(l_j-m)*N_n_j(1)*r(1)/(2*l_j-1) - m*N_l_j(1)*q(1)

         else

            c_t(1) = -(n_j+2)*(n_j+m+1)*N_l_j(1)*q(1)/(2*n_j+3) + m*N_n_j(1)*r(1)
            c_p(1) = (l_j+2)*(l_j+m+1)*N_n_j(2)*r(2)/(2*l_j+3) - (l_j-1)*(l_j-m)*N_n_j(1)*r(1)/(2*l_j-1) - &
                     m*N_l_j(1)*q(1)

            l_j = l_j + 2
            n_j = n_j + 2

            do j = 2,n-1
               c_t(j) = (n_j-1)*(n_j-m)*N_l_j(j-1)*q(j-1)/(2*n_j-1) - (n_j+2)*(n_j+m+1)*N_l_j(j)*q(j)/(2*n_j+3) + &
                        m*N_n_j(j)*r(j)
               c_p(j) = (l_j+2)*(l_j+m+1)*N_n_j(j+1)*r(j+1)/(2*l_j+3) - (l_j-1)*(l_j-m)*N_n_j(j)*r(j)/(2*l_j-1) - &
                        m*N_l_j(j)*q(j)
               l_j = l_j + 2
               n_j = n_j + 2
            end do

            c_t(n) = (n_j-1)*(n_j-m)*N_l_j(n-1)*q(n-1)/(2*n_j-1) - (n_j+2)*(n_j+m+1)*N_l_j(n)*q(n)/(2*n_j+3) + &
                     m*N_n_j(n)*r(n)
            c_p(n) = -(l_j-1)*(l_j-m)*N_n_j(n)*r(n)/(2*l_j-1) - m*N_l_j(n)*q(n)

         endif

       end associate

    end if

    ! Adjust signs by considering the behavior at the pole

    if (parity) then
       P_lm = P_lm_series(2*n, m, 0._WP)
       if (DOT_PRODUCT(P_lm(1::2), c_r) < 0._WP) then
          c_r = -c_r
          c_t = -c_t
          c_p = -c_p
       endif
    else
       dP_lm = dP_lm_series(2*n, m, 0._WP)
       if (DOT_PRODUCT(dP_lm(2::2), c_r) < 0._WP) then
          c_r = -c_r
          c_t = -c_t
          c_p = -c_p
       endif
    endif

    ! Finish

    return

  end subroutine eval_coeffs_

!****

  subroutine eval_eigvecs_ (nu, m, k, lambda, b, q, r)

    real(WP), intent(in)  :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    real(WP), intent(in)  :: lambda
    real(WP), intent(in)  :: b(:)
    real(WP), intent(out) :: q(:)
    real(WP), intent(out) :: r(:)

    integer :: n
    logical :: parity
    integer :: j
    integer :: l_j
    integer :: n_j

    $CHECK_BOUNDS(SIZE(q),SIZE(b))
    $CHECK_BOUNDS(SIZE(r),SIZE(b))

    ! Calculate the q and r eigenvectors associated with the b
    ! eigenvector (see Townsend 1997)

    ! q eigenvector (poloidal)

    n = SIZE(b)

    parity = MOD(k, 2) == 0

    if (parity .AND. m == 0) then

       ! Special casing for parity=.TRUE., m=0 situations, since q(1)
       ! is not properly defined; its actual value doesn't matter,
       ! since it has no physical significance, so we set it to zero

       q(1) = 0._WP

       if(k == 0) then
          
          q(2:) = 0._WP

       else

          l_j = 2

          do j = 2, n
             q(j) = b(j)*lambda/(l_j*(l_j+1))
             l_j = l_j + 2
          end do

       endif

    else

       if (parity) then
          l_j = ABS(m)
       else
          l_j = ABS(m) + 1
       endif
       
       do j = 1, n
          q(j) = b(j)*lambda/(l_j*(l_j+1))
          l_j = l_j + 2
       end do

    endif

    ! r eigenvector (toroidal)

    if (parity) then

       if(m == 0 .AND. k == 0) then

          r = 0._WP

       else

          n_j = ABS(m) + 1

          do j = 1, n-1
             r(j) = nu*(J_lm_(n_j, m)*b(j)/n_j**2 + J_lm_(n_j+1, m)*b(j+1)/(n_j+1)**2)*lambda/(1._WP + m*nu/(n_j*(n_j+1)))
             n_j = n_j + 2
          enddo

          r(n) = nu*(J_lm_(n_j, m)*b(j)/n_j**2)*lambda/(1._WP + m*nu/(n_j*(n_j+1)))

       end if

    else

       n_j = ABS(m)

       if (m == 0) then

          ! Special casing for parity=.FALSE., m=0 situations, since
          ! r(1) is not properly defined; its actual value doesn't
          ! matter, since it has no physical significance, so we set
          ! it to zero

          r(1) = 0._WP

       else

          r(1) = nu*(J_lm_(n_j+1, m)*b(1)/(n_j+1)**2)*lambda/(1._WP + m*nu/(n_j*(n_j+1)))

       endif

       n_j = n_j + 2
       
       do j = 2,n
          r(j) = nu*(J_lm_(n_j+1, m)*b(j)/(n_j+1)**2 + J_lm_(n_j, m)*b(j-1)/n_j**2)*lambda/(1._WP + m*nu/(n_j*(n_j+1)))
          n_j = n_j + 2
       end do
       
    endif

    ! Finish

    return

  end subroutine eval_eigvecs_

!****

  subroutine eigen_W_n_ (nu, m, k, n, lambda, b)

    real(WP), intent(in)                         :: nu
    integer, intent(in)                          :: m
    integer, intent(in)                          :: k
    integer, intent(in)                          :: n
    real(WP), intent(out)                        :: lambda
    real(WP), allocatable, optional, intent(out) :: b(:)

    $if ($IEEE)
    logical  :: h_over
    logical  :: h_zero
    logical  :: h_inval
    $endif
    logical  :: parity
    real(WP) :: lambda_A(n)
    real(WP) :: b_A(n,n)
    integer  :: i

    ! Calculate the W matrix eigenvalue and (possibly) eigenvector for
    ! the supplied mode parameters and the indicated matrix dimension

    $if ($IEEE)
    call save_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Calculate the full eigenvalue/eigenvector spectrum

    parity = MOD(k, 2) == 0

    if (PRESENT(b)) then
       call eigen_A_(nu, m, parity, n, lambda_A, b_A)
    else
       call eigen_A_(nu, m, parity, n, lambda_A)
    endif

    ! Determine the index of the required eigenpair

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

    ! Pick out the eigenpair

    lambda = 1._WP/lambda_A(i)
    if (PRESENT(b)) b = b_A(:,i)

    $if ($IEEE)
    call restore_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Finish

    return

  end subroutine eigen_W_n_

!****

  subroutine eigen_W_tol_ (nu, m, k, tol, lambda, b)

    real(WP), intent(in)                         :: nu
    integer, intent(in)                          :: m
    integer, intent(in)                          :: k
    real(WP), intent(in)                         :: tol
    real(WP), intent(out)                        :: lambda
    real(WP), allocatable, optional, intent(out) :: b(:)

    real(WP), parameter :: EPS = EPSILON(0._WP)
    
    logical               :: parity
    integer               :: n
    real(WP)              :: lambda_cmp
    real(WP), allocatable :: b_cmp(:)

    ! Calculate the W matrix eigenvalue and (possibly) eigenvector for
    ! the supplied mode parameters and the indicated tolerance

    ! Initialize n to the smallest possible value

    parity = MODULO(k, 2) == 0

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

    ! Now increase n until the desired tolerance is reached

!    allocate(b(n))

    call eigen_W_n_(nu, m, k, n, lambda, b)

    if (PRESENT(b)) b = b/MAXVAL(ABS(b))

    n_loop : do

       ! if (PRESENT(b)) then

       !    allocate(b_cmp(2*n))

       !    call eigen_W_n_(nu, m, k, 2*n, lambda_cmp, b_cmp)
       !    b_cmp = b_cmp/MAXVAL(ABS(b_cmp))

       !    print *,'Trying:',nu,n,COUNT(ABS(b - b_cmp(:n)) <= (4._WP*EPS+tol)*ABS(b_cmp(:n)))

       !    if(ALL(ABS(b - b_cmp(:n)) <= (4._WP*EPS+tol)*ABS(b_cmp(:n)))) exit n_loop

       !    lambda = lambda_cmp
       !    b = b_cmp

       !    deallocate(b_cmp)

       ! else

          call eigen_W_n_(nu, m, k, 2*n, lambda_cmp)

          if(ABS(lambda - lambda_cmp) <= (4._WP*EPS+tol)*ABS(lambda_cmp)) exit n_loop

          lambda = lambda_cmp

       ! endif

       n = 2*n

    end do n_loop

    ! Finish

    return

  end subroutine eigen_W_tol_

!****

  subroutine eigen_A_ (nu, m, parity, n, lambda, b)

    real(WP), intent(in)            :: nu
    integer, intent(in)             :: m
    logical, intent(in)             :: parity
    integer, intent(in)             :: n
    real(WP), intent(out)           :: lambda(:)
    real(WP), optional, intent(out) :: b(:,:)

    $if ($IEEE)
    logical  :: h_over
    logical  :: h_zero
    logical  :: h_inval
    $endif
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
    real(WP) :: abstol
    integer  :: j
    integer  :: j_inf

    $CHECK_BOUNDS(SIZE(lambda),n)

    if (PRESENT(b)) then
       $CHECK_BOUNDS(SIZE(b, 1),n)
       $CHECK_BOUNDS(SIZE(b, 2),n)
    endif

    ! Calculate the A (alpeh) matrix eigenvalues and (possibly)
    ! eigenvectors for the supplied mode parameters

    $if ($IEEE)
    call save_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Assemble the A-matrix

    call assemble_A_(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, &
                     A_D, A_E)
    
    A_E(n) = 0._WP

    ! Calculate the eigenvalues and eigenvectors, taking care to
    ! handle special cases where one or more components of A is
    ! infinite

    abstol = LA_LAMCH(1.0_WP, 'Safe minimum')

    if (A_5(1) == 0._WP) then

       ! A(1,1) infinite, causing a +inf eigenvalue

       if (PRESENT(b)) then

          call LA_STEVR(A_D(2:), A_E(2:), lambda(:n-1), Z=b(:n-1,:n-1))

          b(n,:n-1) = 0._WP
          b(:n-1,n) = 0._WP

          b(n,n) = 1._WP

       else

          call LA_STEVR(A_D(2:), A_E(2:), lambda(:n-1))

       endif

       $if ($IEEE)
       lambda(n) = IEEE_VALUE(0._WP, IEEE_POSITIVE_INF)
       $else
       lambda(n) = HUGE(0._WP)
       $endif

    elseif (A_2d(1) == 0._WP .AND. A_2n(1) /= 0._WP) then

       ! A(1,1) infinite, causing a -inf eigenvalue

       if (PRESENT(b)) then

          call LA_STEVR(A_D(2:), A_E(2:), lambda(2:), Z=b(2:,2:))

          b(1,2:) = 0._WP
          b(2:,1) = 0._WP

          b(1,1) = 1._WP

       else

          call LA_STEVR(A_D(2:), A_E(2:), lambda(2:))

       endif

       $if ($IEEE)
       lambda(1) = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)
       $else
       lambda(1) = -HUGE(0._WP)
       $endif

    elseif (A_3d(n) == 0._WP .AND. A_3n(n) /= 0._WP) then

       ! A(n,n) infinite, causing a -inf eigenvalue

       if (PRESENT(b)) then

          call LA_STEVR(A_D(:n-1), A_E(:n-1), lambda(2:), Z=b(2:,2:))

          b(1,2:) = 0._WP
          b(2:,1) = 0._WP

          b(1,1) = 1._WP

       else

          call LA_STEVR(A_D(:n-1), A_E(:n-1), lambda(2:))

       endif

       $if ($IEEE)
       lambda(1) = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)
       $else
       lambda(1) = -HUGE(0._WP)
       $endif

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

          $ABORT(Pure Rossby mode encountered)

       else

          call LA_STEVR(A_D, A_E, lambda, Z=b)

       endif

    endif

    $if ($IEEE)
    call restore_ieee_state_(h_over, h_zero, h_inval)
    $endif

    ! Finish

    return

  end subroutine eigen_A_

!****

  subroutine eval_A_ (nu, m, n, parity, A_D, A_E)

    real(WP), intent(in)  :: nu
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    logical, intent(in)    :: parity
    real(WP), intent(out) :: A_D(:)
    real(WP), intent(out) :: A_E(:)

    real(WP) :: A_1(n)
    real(WP) :: A_2n(n)
    real(WP) :: A_2d(n)
    real(WP) :: A_3n(n)
    real(WP) :: A_3d(n)
    real(WP) :: A_4n(n)
    real(WP) :: A_4d(n)
    real(WP) :: A_5(n)
    integer   :: j

    $CHECK_BOUNDS(SIZE(A_D),n)
    $CHECK_BOUNDS(SIZE(A_E),n)

    ! Evaluate the A matrix

    call assemble_A_(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, A_D, A_E)

    ! Finish

    return

  end subroutine eval_A_

!****

  subroutine assemble_A_ (nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, A_D, A_E)

    real(WP), intent(in)  :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: n
    logical, intent(in)   :: parity
    real(WP), intent(out) :: A_1(:)
    real(WP), intent(out) :: A_2n(:)
    real(WP), intent(out) :: A_2d(:)
    real(WP), intent(out) :: A_3n(:)
    real(WP), intent(out) :: A_3d(:)
    real(WP), intent(out) :: A_4n(:)
    real(WP), intent(out) :: A_4d(:)
    real(WP), intent(out) :: A_5(:)
    real(WP), intent(out) :: A_D(:)
    real(WP), intent(out) :: A_E(:)

    integer  :: j
    integer  :: l_j
    real(WP) :: A_2
    real(WP) :: A_3
    real(WP) :: A_4

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

    if (parity) then
       l_j = ABS(m)
    else
       l_j = ABS(m) + 1
    endif
        
    do j = 1, n

       if (l_j == ABS(m)) then
          A_1(j) = 1._WP + SIGN(1, m)*nu/REAL(l_j+1, WP)
          A_2n(j) = 0._WP
       else
          A_1(j) = 1._WP + m*nu/(REAL(l_j, WP)*REAL(l_j+1, WP))
          A_2n(j) = -nu**2*REAL(l_j-1, WP)**2*REAL(l_j+1, WP)*J_lm_(l_j, m)**2/REAL(l_j, WP)
       endif

       A_2d(j) = REAL(l_j-1, WP)*REAL(l_j, WP) + m*nu
       
       A_3n(j) = -nu**2*REAL(l_j, WP)*REAL(l_j+2, WP)**2*J_lm_(l_j+1, m)**2/REAL(l_j+1, WP)
       A_3d(j) = REAL(l_j+1, WP)*REAL(l_j+2, WP) + m*nu

       A_4n(j) = -nu**2*J_lm_(l_j+1, m)*J_lm_(l_j+2, m)
       A_4d(j) = A_3d(j)

       A_5(j) = REAL(l_j, WP)*REAL(l_j+1, WP)

       l_j = l_j + 2

    end do

    ! Now assemble the matrix

    do j = 1, n-1

       if (A_2n(j) /= 0._WP) then
          A_2 = A_2n(j)/A_2d(j)
       else
          A_2 = 0._WP
       endif

       if (A_3n(j) /= 0._WP) then
          A_3 = A_3n(j)/A_3d(j)
       else
          A_3 = 0._WP
       endif

       if (A_4n(j) /= 0._WP) then
          A_4 = A_4n(j)/A_4d(j)
       else
          A_4 = 0._WP
       endif

       if (A_5(j) /= 0._WP) then
          A_D(j) = (A_1(j) + A_2 + A_3)/A_5(j)
       endif
       A_E(j) = A_4

    end do

    if (A_2n(n) /= 0._WP) then
       A_2 = A_2n(n)/A_2d(n)
    else
       A_2 = 0._WP
    endif

    if (A_3n(n) /= 0._WP .AND. parity) then
       A_3 = A_3n(n)/A_3d(n)
    else
       A_3 = 0._WP
    endif

    if (A_5(n) /= 0._WP) then
       A_D(n) = (A_1(n) + A_2 + A_3)/A_5(n)
    endif

    ! Finish

    return

  end subroutine assemble_A_

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

    if(-m*nu > 0._WP) then

       if(parity) then
          l_j = ABS(m)
       else
          l_j = ABS(m) + 1
       endif
       
       j_loop : do j = 1,n
          if (-m*nu > REAL(l_j+1, WP)*REAL(l_j+2, WP)) Xi_ = Xi_ + 1
          l_j = l_j + 2
       enddo j_loop
       
       if(.NOT. parity .AND. -m*nu > REAL(ABS(m), WP)*REAL(ABS(m)+1, WP)) Xi_ = Xi_ + 1
       
    endif

    ! Finish

    return

  end function Xi_

!****

  function J_lm_ (l, m)

    integer, intent(in) :: l
    integer, intent(in) :: m
    real(WP)            :: J_lm_

    ! Calculate the J_lm function, being the integral of Y_lm Y_{l+1}m
    ! cos(theta) over all solid angles

    if(ABS(m) >= l) then

       J_lm_ = 0._WP

    else

       J_lm_ = SQRT(REAL(l-m, WP)*REAL(l+m, WP)/(REAL(2*l+1, WP)*REAL(2*l-1, WP)))

    endif

    ! Finish

    return
    
  end function J_lm_

!****

  function P_lm_series (n, m, x) result (P_lm)

    integer, intent(in)  :: n
    integer, intent(in)  :: m
    real(WP), intent(in) :: x
    real(WP)             :: P_lm(n)

    integer  :: am
    real(WP) :: sqrt_term
    real(WP) :: fact_term
    integer  :: i
    integer  :: l

    ! Calculate an n-term series of associated Legendre functions
    ! P_lm(x), with l=|m|,|m|+1,...,|m|+n-1

    if (n > 0) then

       ! First calculate the P_mm term

       P_lm(1) = 1._WP
       am = ABS(m)

       if (m > 0) then

          sqrt_term = -SQRT((1._WP-x)*(1._WP+x))
          fact_term = 1._WP

          do i = 1,am
             P_lm(1) = P_lm(1)*fact_term*sqrt_term
             fact_term = fact_term + 2._WP
          end do
       
       elseif (m < 0) then

          sqrt_term = SQRT((1._WP-x)*(1._WP+x))
          fact_term = 2._WP

          do i = 1,am
             P_lm(1) = P_lm(1)*sqrt_term/fact_term
             fact_term = fact_term + 2.
          end do

       endif

       ! Now calculate the other terms using recursion

       if (n > 1) then

          P_lm(2) = x*(2*am+1)*P_lm(1)/(am+1-m)

          l = am + 2

          do i = 3,n
             P_lm(i) = (x*(2*l-1)*P_lm(i-1) - (l+m-1)*P_lm(i-2))/(l-m)
             l = l + 1
          end do

       endif

    end if

    ! Finish

    return

  end function P_lm_series

!****

  function dP_lm_series (n, m, x) result (dP_lm)

    integer, intent(in)  :: n
    integer, intent(in)  :: m
    real(WP), intent(in) :: x
    real(WP)             :: dP_lm(n)

    integer  :: am
    real(WP) :: P_lm_sqrt(n)
    real(WP) :: sqrt_term
    real(WP) :: fact_term
    integer  :: i
    integer  :: l

    ! Calculate an n-term series of associated Legendre function derivatives
    ! -sqrt(1-x^2) dP_lm(x)/dx = dP_lm(cos theta)/dtheta, with l=|m|,|m|+1,...,|m|+n-1

    if(n > 0) then

       if (m /= 0) then

          ! This algorithm is based on calculating the functions
          ! P_lm/sqrt(1-x^2), and then using the formula for (x^2-1)
          ! dP_lm/dx

          ! First calculate the P_mm term

          am = ABS(m)

          if (m > 0) then

             P_lm_sqrt(1) = 1._WP

             sqrt_term = -SQRT((1._WP-x)*(1._WP+x))
             fact_term = 3._WP

             do i = 2,am
                P_lm_sqrt(1) = P_lm_sqrt(1)*fact_term*sqrt_term
                fact_term = fact_term + 2._WP
             end do
       
          elseif (m < 0) then

             P_lm_sqrt(1) = -0.5_WP

             sqrt_term = SQRT((1._WP-x)*(1._WP+x))
             fact_term = 4._WP

             do i = 2,am
                P_lm_sqrt(1) = P_lm_sqrt(1)*sqrt_term/fact_term
                fact_term = fact_term + 2._WP
             end do

          endif
          
          ! Now calculate the other terms using recursion

          if (n > 1) then

             P_lm_sqrt(2) = x*(2*am+1)*P_lm_sqrt(1)/(am+1-m)

             l = am + 2

             do i = 3,n
                P_lm_sqrt(i) = (x*(2*l-1)*P_lm_sqrt(i-1) - (l+m-1)*P_lm_sqrt(i-2))/(l-m)
                l = l + 1
             end do

          endif

          ! Apply the formula

          dP_lm(1) = am*x*P_lm_sqrt(1)

          l = am + 1

          do i = 2, n
             dP_lm(i) = l*x*P_lm_sqrt(i) - (l+m)*P_lm_sqrt(i-1)
             l = l + 1
          end do

       else

          ! m=0 is a special case

          dP_lm(1) = 0.
          dP_lm(2:) = -P_lm_series(n-1, 1, x)

       endif

    endif

    ! Finish

    return

  end function dP_lm_series

!****

  function N_lm_series (n, m) result (N_lm)

    integer, intent(in) :: n
    integer, intent(in) :: m
    real(WP)            :: N_lm(n)

    integer  :: am
    real(WP) :: fact_term
    integer  :: i
    integer  :: l

    ! Calculate an n-term series of spherical harmonic normalizing
    ! factors N_lm, with l=|m|,|m|+1,...,|m|+n-1

    if (n > 0) then

       am = ABS(m)

       fact_term = 1._WP
    
       if (m >= 0) then

          do i = 1,2*am
             fact_term = fact_term/i
          end do

       else

          do i = 1,2*am
             fact_term = fact_term*i
          end do

       endif

       N_lm(1) = fact_term*(2*am+1)

       l = am + 1

       if (m >= 0) then

          do i = 2,n
             fact_term = (l-am)*fact_term/(l+am)
             N_lm(i) = fact_term*(2*l+1)
             l = l + 1
          end do

       else

          do i = 2,n
             fact_term = (l+am)*fact_term/(l-am)
             N_lm(i) = fact_term*(2*l+1)
             l = l + 1
          end do

       endif

       N_lm = SQRT(0.5_WP)*SQRT(N_lm)

    endif

    ! Finish

    return

  end function N_lm_series

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

!****

  function lambda_n_ (nu, m, k, n) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    integer, intent(in)  :: n
    real(WP)             :: lambda
    
    ! Evaluate the Hough eigenvalue lambda for the supplied mode
    ! parameters and the indicated matrix dimension

    call eigen_W_n_(nu, m, k, n, lambda)

    ! Finish

    return

  end function lambda_n_

!****

  function lambda_tol_ (nu, m, k, tol) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: tol
    real(WP)             :: lambda
    
    ! Evaluate the Hough eigenvalue lambda for the supplied mode
    ! parameters and the indicated tolerance

    call eigen_W_tol_(nu, m, k, tol, lambda)

    ! Finish

    return

  end function lambda_tol_

!****

  function lambda_asymp (nu, m, k) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda

    integer :: s
    
    ! Evaluate the (m, k) Hough eigenvalue lambda using the asymptotic
    ! expressions by Townsend (2003, MNRAS, 340,1020)

    if (m*nu >= 0) then

       $ASSERT(k >= 0,Invalid k)

       if (k > 0) then

          s = k - 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))


       else

          lambda = m**2*(2._WP*m*nu)/(2._WP*m*nu - 1._WP)

       endif

    else

       if (k >= -1) then

          s = k + 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))

       else

          s = -k -1

          lambda = (m*nu - m**2)**2/(nu**2*(2*s+1)**2)

       endif

    endif

    ! Finish

    return

  end function lambda_asymp

end module astro_hough
