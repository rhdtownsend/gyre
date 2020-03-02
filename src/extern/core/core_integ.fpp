  ! Module : core_integ
  ! purpose : numerical integration olio

$include 'core.inc'

module core_integ

  ! Uses

  use core_kinds
  use core_constants
  
  use ISO_FORTRAN_ENV

  ! no implicit typing
  
  implicit none

  ! Access specifiers

  private

  public :: simps2d_wgt
  public :: simps1d_wgt
  public :: gl_weights

contains
  
 !****

  subroutine simps2d_wgt(wgt, i, j, ni, nj)

    real(WP), intent(out)    :: wgt
    integer, intent(in)      :: i
    integer, intent(in)      :: j
    integer, intent(in)      :: ni
    integer, intent(in)      :: nj

    real(WP)  :: wgt1
    real(WP)  :: wgt2

    call simps1d_wgt(wgt1, i, ni)
    call simps1d_wgt(wgt2, j, nj)
    
    wgt = wgt1 * wgt2

    return
    
  end subroutine simps2d_wgt
  
  !****

  subroutine simps1d_wgt(wgt, i, n)
    real(WP), intent(out)    :: wgt
    integer, intent(in)      :: i
    integer, intent(in)      :: n

    $ASSERT(i <= n, wrong index in simps1d)
    $ASSERT(MOD(n,2) == 1, n must be odd)
    
    if(i == 1 .OR. i == n) then
       wgt = 1._WP
    else if(MOD(i,2) == 1) then
       wgt = 2._WP
    else
       wgt = 4._WP
    end if
        
    return
    
  end subroutine simps1d_wgt

  !****

  ! from Numerical Receipes, Press et al
  ! calculate abyssica (x) and weights (w) for Gauss-Legendre quadrature
  ! over integral from a to b
  
  subroutine gl_weights(a, b, x, w)
    real(WP), intent(in) :: a
    real(WP), intent(in) :: b
    real(WP), intent(inout) :: x(:)
    real(WP), intent(inout) :: w(:)

    integer :: nx
    integer :: m
    integer :: i
    integer :: j
    real(WP) :: z1
    real(WP) :: z
    real(WP) :: xm
    real(WP) :: xl
    real(WP) :: pp
    real(WP) :: p3
    real(WP) :: p2
    real(WP) :: p1 
    
    nx = SIZE(x)

    $ASSERT(nx == SIZE(w), wrong size of input x or w)

    m = (nx+1)/2

    xm = 0.5_WP*(b+a)
    xl = 0.5_WP*(b-a)

    do i=1,m
       z = COS(PI*((i-1)+0.75_WP)/(nx+0.5_WP))

       do
          p1 = 1._WP
          p2 = 0._WP

          do j=1,nx
             p3 = p2
             p2 = p1
             p1 = ((2._WP*(j-1)+1._WP)*z*p2-(j-1)*p3)/(j)
          end do

          pp = nx*(z*p1-p2)/(z**2-1._WP)
          z1 = z
          z = z1-p1/pp

          if(ABS(z-z1) < 10._WP*EPSILON(0._WP)) then
             exit
          end if
          
       end do

       x(i) = xm - xl*z
       x(nx-i+1) = xm + xl*z
       w(i) = 2._WP*xl/((1._WP-z**2)*pp**2)
       w(nx-i+1) = w(i)
    end do

    return

  end subroutine gl_weights
  
end module core_integ

