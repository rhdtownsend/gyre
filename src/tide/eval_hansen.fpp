program eval_hansen

  ! Uses

  use core_kinds
  use core_system

  use gyre_math
  use gyre_tide_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP) :: e
  integer  :: n
  integer  :: m
  integer  :: k
  real(WP) :: X
  real(WP) :: X_hat
  real(WP) :: X_hat_ck
  real(WP) :: X_tilde
  real(WP) :: X_tilde_ck

  ! Get arguments

  if (n_arg() /= 4) stop '** syntax: eval_hansen e n m k'

  call get_arg(1, e)
  call get_arg(2, n)
  call get_arg(3, m)
  call get_arg(4, k)

  ! Initialize

  call init_math()

  ! Print out the hansen coefficient

  X = hansen_X(e, n, m, k)

  X_hat = hansen_X_hat(e, n, m, k)
  X_tilde = hansen_X_tilde(e, n, m, k)

  X_hat_ck = 0.5_WP*(hansen_X(e, n, m-1, k) + hansen_X(e, n, m+1, k))
  X_tilde_ck = 0.5_WP*(hansen_X(e, n, m-1, k) - hansen_X(e, n, m+1, k))

!  print *, X, X_hat, X_hat_ck, X_tilde, X_tilde_ck
  print *, X

end program eval_hansen
  
