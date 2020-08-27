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

  integer  :: n_eval

  integer  :: i
  real(WP) :: e
  real(WP) :: r
  integer  :: n
  integer  :: m
  integer  :: k
  real(WP) :: X

  integer :: c_beg, c_end, c_rate, c_hansen

  ! Get arguments

  if (n_arg() /= 1) stop '** syntax: time_eval_hansen n_eval'

  call get_arg(1, n_eval)

  ! Initialize

  call init_math()

  ! Loop over evaluations

  do i = 1, n_eval

     ! Generate randomized numbers

     call RANDOM_NUMBER(e)

     call RANDOM_NUMBER(r)
     n = -10 + 20*r

     call RANDOM_NUMBER(r)
     m = -10 + 20*r

     call RANDOM_NUMBER(r)
     k = -10 + 20*r

!     print *,e, n, m, k

     call SYSTEM_CLOCK(c_beg, c_rate)

     X = hansen_X(e, n, m, k)

     call SYSTEM_CLOCK(c_end)

     c_hansen = c_hansen + (c_end - c_beg)

  end do

  print *,'Timing:',REAL(c_hansen)/c_rate

end program eval_hansen
  
