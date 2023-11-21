$include 'core.inc'

program test_interp

  use kinds_m
  
  use interp_m

  real(WP), parameter :: x_a = -3._WP
  real(WP), parameter :: x_b = 3._WP
  integer, parameter  :: n = 20
  integer, parameter  :: m = 100

  type(r_interp_t) :: in
  real(WP)         :: x(n)
  real(WP)         :: y(n)
  real(WP)         :: dy_dx_a
  real(WP)         :: dy_dx_b
  real(WP)         :: xx

  ! Set up x and y

  do i = 1, n
     x(i) = ((n-i)*x_a + (i-1)*x_b)/(n-1)
     y(i) = 3._WP*x(i)**3 - x(i) + 4._WP
  end do

  dy_dx_a = 9._WP*x(1)**2 - 1._WP
  dy_dx_b = 9._WP*x(n)**2 - 1._WP

  ! Set up the interpolant

  in = r_interp_t_eval_derivs_(x, y, 'SPLINE', df_dx_a=dy_dx_a, df_dx_b=dy_dx_b)

  ! Now calculate

  do i = 1, m
     xx = ((m-i)*x_a + (i-1)*x_b)/(m-1)
     print *, xx, in%f(xx), in%df_dx(xx, 0), in%df_dx(xx), in%df_dx(xx, 2), in%df_dx(xx, 3), in%int_f(xx)
  end do

  ! Finish

end program test_interp
