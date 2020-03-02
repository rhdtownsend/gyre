! Program  : astro_hough
! Purpose  : calculation of Hough eigenvalues & eigenfunctions

$include 'core.inc'

module astro_hough

  ! Uses

  use core_kinds
  use core_linalg
  use core_order

  use ISO_FORTRAN_ENV

  use IEEE_ARITHMETIC
  
  ! No implicit typing

  implicit none

  ! Interfaces

  interface lambda_fix
     module procedure lambda_fix_r_
  end interface lambda_fix

  interface lambda_var
     module procedure lambda_var_r_
  end interface lambda_var

  ! Access specifiers

  private

  public :: lambda_asymp
  public :: lambda_fix
  public :: lambda_var

  ! Procedures

contains

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

       else

          lambda = m**2*(2._WP*m*nu)/(2._WP*m*nu - 1._WP)

       endif

    else

       if (k >= -1) then

          s = k + 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

       else

          s = -k -1

          lambda = (m*nu - m**2)**2/(nu**2*(2*s+1)**2)

       endif

    endif

    ! Finish

    return

  end function lambda_asymp

!****

  function lambda_fix_r_ (nu, m, k, n) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    integer, intent(in)  :: n
    real(WP)             :: lambda

    logical  :: h_over
    logical  :: h_zero
    logical  :: h_inval
    logical  :: parity
    real(WP) :: y(n)
    integer  :: i

    ! Evaluate the (m,k) Hough eigenvalue lambda from an A-matrix with
    ! fixed dimension n

    call IEEE_GET_HALTING_MODE(IEEE_OVERFLOW, h_over)
    call IEEE_GET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, h_zero)
    call IEEE_GET_HALTING_MODE(IEEE_INVALID, h_inval)

    call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE.)
    call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE.)
    call IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)

    ! Calculate the eigenvalues of A

    parity = MOD(k, 2) == 0

    call eval_y_r_(nu, m, n, parity, y)

    ! Determine the index of the required eigenvalue

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

    ! Calculate the Hough eigenvalue

    lambda = 1._WP/y(i)

    ! Finish

    call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, h_over)
    call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, h_zero)
    call IEEE_SET_HALTING_MODE(IEEE_INVALID, h_inval)

    return

  end function lambda_fix_r_

!****

  function lambda_var_r_ (nu, m, k, lambda_tol) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: lambda_tol
    real(WP)             :: lambda
    
    real(WP), parameter :: EPS = EPSILON(0._WP)

    logical  :: parity
    integer  :: n
    real(WP) :: lambda_cmp

    ! Evaluate the (m,k) Hough eigenvalue lambda from an A-matrix with
    ! variable dimension, adjusted to give tolerance lambda_tol

    ! Initialize n to the smallest possible value

    parity = MODULO(k, 2) == 0

    if(k >= 0) then

       if(parity) then
          n = k/2+1
       else
          n = (k+1)/2+1
       endif

       if (m*nu < 0._WP) then
          n = n + Xi_(nu, m, n, parity)
       endif

    else

       if(parity) then
          n = ABS(k)/2
       else
          n = ABS(k-1)/2
       endif

    endif

    ! Now increase n until the desired tolerance is reached

    lambda = lambda_fix(nu, m, k, n)

    n_loop : do

       lambda_cmp = lambda_fix(nu, m, k, 2*n)

       if(ABS(lambda - lambda_cmp) <= (4._WP*EPS+lambda_tol)*ABS(lambda_cmp)) exit n_loop

       n = 2*n

       lambda = lambda_cmp

    end do n_loop

    ! Finish

    return

  end function lambda_var_r_

!****

  $define $EVAL_A $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine eval_A_${INFIX}_ (nu, m, n, parity, A_D, A_E)

    $TYPE(WP), intent(in)  :: nu
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    logical, intent(in)    :: parity
    $TYPE(WP), intent(out) :: A_D(:)
    $TYPE(WP), intent(out) :: A_E(:)

    $TYPE(WP) :: A_1(n)
    $TYPE(WP) :: A_2n(n)
    $TYPE(WP) :: A_2d(n)
    $TYPE(WP) :: A_3n(n)
    $TYPE(WP) :: A_3d(n)
    $TYPE(WP) :: A_4n(n)
    $TYPE(WP) :: A_4d(n)
    $TYPE(WP) :: A_5(n)
    integer   :: j
    $TYPE(WP) :: A_2
    $TYPE(WP) :: A_3
    $TYPE(WP) :: A_4

    $CHECK_BOUNDS(SIZE(A_D),n)
    $CHECK_BOUNDS(SIZE(A_E),n-1)

    ! Evaluate the A matrix

    call assemble_A(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, A_D, A_E)

    ! Finish

    return

  end subroutine eval_A_${INFIX}_

  $endsub

  $EVAL_A(r,real)
  $EVAL_A(c,complex)

!****

  subroutine eval_y_r_ (nu, m, n, parity, y)

    real(WP), intent(in)  :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: n
    logical, intent(in)   :: parity
    real(WP), intent(out) :: y(:)

    real(WP) :: A_1(n)
    real(WP) :: A_2n(n)
    real(WP) :: A_2d(n)
    real(WP) :: A_3n(n)
    real(WP) :: A_3d(n)
    real(WP) :: A_4n(n)
    real(WP) :: A_4d(n)
    real(WP) :: A_5(n)
    real(WP) :: A_D(n)
    real(WP) :: A_E(n-1)
    real(WP) :: abstol
    integer  :: n_y
    real(WP) :: Z(1,1)
    integer  :: isuppz(n)
    real(WP) :: work(20*n)
    integer  :: iwork(10*n)
    integer  :: info
    logical  :: inf_block
    integer  :: j
    
    $CHECK_BOUNDS(SIZE(y),n)

    ! Evaluate the eigenvalues y of the A matrix

    ! Calculate A

    call assemble_A_r_(nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, &
                       A_D, A_E)

    ! Calculate the eigenvalues, applying heuristics to handle
    ! infinite components

    abstol = XLAMCH(0._WP, 'S')

    if (A_5(1) == 0._WP) then

       ! A(1,1) infinite

       call DSTEVR('N', 'A', n-1, A_D(2:), A_E(2:), 0._WP, 0._WP, 0, 0, abstol, n_y, &
                   y(:n-1), Z, n, isuppz, work, 20*n, iwork, 10*n, info)

       y(n) = IEEE_VALUE(0._WP, IEEE_POSITIVE_INF)

    elseif (A_2d(1) == 0._WP .AND. A_2n(1) /= 0._WP) then

       ! A(1,1) infinite

       call DSTEVR('N', 'A', n-1, A_D(2:), A_E(2:), 0._WP, 0._WP, 0, 0, abstol, n_y, &
                   y(2:), Z, n, isuppz, work, 20*n, iwork, 10*n, info)

       y(1) = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)

    elseif (A_3d(n) == 0._WP .AND. A_3n(n) /= 0._WP) then

       ! A(n,n) infinite

       call DSTEVR('N', 'A', n-1, A_D(:n-1), A_E(:n-1), 0._WP, 0._WP, 0, 0, abstol, n_y, &
                   y(2:), Z, n, isuppz, work, 20*n, iwork, 10*n, info)

       y(1) = IEEE_VALUE(0._WP, IEEE_NEGATIVE_INF)

    else
       
       ! Scan for a 2x2 block of infinite components

       inf_block = .FALSE.

       scan_loop : do j = 1, n-1
          if (A_3d(j) == 0._WP .AND. A_3n(j) /= 0._WP) then
             inf_block = .TRUE.
             exit scan_loop
          endif
       end do scan_loop

       if (inf_block) then

          ! Rebuild A with a nearby nu (it would be better to split out the bad components...)

          call assemble_A_r_(nu*(1._WP-SQRT(EPSILON(0._WP))), m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, &
                             A_D, A_E)
          
       endif

       call DSTEVR('N', 'A', n, A_D, A_E, 0._WP, 0._WP, 0, 0, abstol, n_y, &
                   y, Z, n, isuppz, work, 20*n, iwork, 10*n, info)

    endif

    ! Finish

    return

  end subroutine eval_y_r_

!****

  $define $ASSEMBLE_A $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine assemble_A_${INFIX}_ (nu, m, n, parity, A_1, A_2n, A_2d, A_3n, A_3d, A_4n, A_4d, A_5, A_D, A_E)

    real(WP), intent(in)   :: nu
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    logical, intent(in)    :: parity
    $TYPE(WP), intent(out) :: A_1(:)
    $TYPE(WP), intent(out) :: A_2n(:)
    $TYPE(WP), intent(out) :: A_2d(:)
    $TYPE(WP), intent(out) :: A_3n(:)
    $TYPE(WP), intent(out) :: A_3d(:)
    $TYPE(WP), intent(out) :: A_4n(:)
    $TYPE(WP), intent(out) :: A_4d(:)
    $TYPE(WP), intent(out) :: A_5(:)
    $TYPE(WP), intent(out) :: A_D(:)
    $TYPE(WP), intent(out) :: A_E(:)

    integer   :: j
    integer   :: l_j
    $TYPE(WP) :: A_2
    $TYPE(WP) :: A_3
    $TYPE(WP) :: A_4

    $CHECK_BOUNDS(SIZE(A_1),n)
    $CHECK_BOUNDS(SIZE(A_2n),n)
    $CHECK_BOUNDS(SIZE(A_2d),n)
    $CHECK_BOUNDS(SIZE(A_3n),n)
    $CHECK_BOUNDS(SIZE(A_3d),n)
    $CHECK_BOUNDS(SIZE(A_4n),n)
    $CHECK_BOUNDS(SIZE(A_4d),n)
    $CHECK_BOUNDS(SIZE(A_5),n)

    $CHECK_BOUNDS(SIZE(A_D),n)
    $CHECK_BOUNDS(SIZE(A_E),n-1)

    ! Evaluate constituent parts of the A matrix. These are defined so
    ! that:
    !
    ! A(j,j)   = [A_1(j) + A_2n(j)/A_2d(j) + A_3n(j)/A_3d(j)] / A_5(j)
    ! A(j+1,j) = A(j, j+1) = A_4n(j)/A_4d(j)

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

       A_D(j) = (A_1(j) + A_2 + A_3)/A_5(j)
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

    A_D(n) = (A_1(n) + A_2 + A_3)/A_5(n)

    ! Finish

    return

  end subroutine assemble_A_${INFIX}_

  $endsub

  $ASSEMBLE_A(r,real)
  $ASSEMBLE_A(c,complex)

!****

  function Xi_ (nu, m, n, parity)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    logical, intent(in)  :: parity
    integer              :: Xi_

    integer :: l_j
    integer :: j

    ! Calculate the permutation constant Xi (see Townsend, 1997,
    ! MNRAS, 284, 839). Note that this code adopts the sign convention
    ! m*nu > 0 for prograde modes

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
  
end module astro_hough
