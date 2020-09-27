! Program  : build_tar_fit
! Purpose  : evaluate TAR eigenvalue
!
! Copyright 2020 Rich Townsend & The GYRE Team
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

program eval_lambda

  ! Uses

  use core_kinds
  use core_hgroup
  use core_system
  use core_parallel

  use gyre_constants
  use gyre_math
  use gyre_tar_eigen

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                   :: m
  integer                   :: l
  real(WP)                  :: q_min
  real(WP)                  :: q_max
  integer                   :: n_q
  logical                   :: log_q
  logical                   :: rossby
  character(:), allocatable :: filename

  integer               :: k
  real(WP), allocatable :: q(:)
  real(WP), allocatable :: lam(:)
  integer               :: i
  type(hgroup_t)        :: hg

  ! Read parameters

  $ASSERT(n_arg() == 8,Syntax: eval_lambda l m q_min q_max n_q log_q rossby filename)

  call get_arg(1, l)
  call get_arg(2, m)
  call get_arg(3, q_min)
  call get_arg(4, q_max)
  call get_arg(5, n_q)
  call get_arg(6, log_q)
  call get_arg(7, rossby)
  call get_arg(8, filename)

  ! Initialize

  call init_parallel()

  call init_math()

  ! Check arguments & set up k

  $ASSERT(ABS(m) <= l,Invalid m)

  if (rossby) then
     k = -(l - ABS(m) + 1)
  else
     k = l - ABS(m)
  endif

  ! Allocate arrays

  allocate(q(n_q))
  allocate(lam(n_q))

  ! Evaluate them

  if (n_q > 1) then

     !$OMP PARALLEL DO SCHEDULE (GUIDED)
     do i = 1, n_q

        if (log_q) then
           q(i) = 10**((LOG10(q_min)*(n_q-i) + LOG10(q_max)*(i-1))/(n_q-1))
        else
           q(i) = (q_min*(n_q-i) + q_max*(i-1))/(n_q-1)
        endif

        lam(i) = lambda(q(i), m, k)

     end do

  elseif (n_q == 1) then

     $ASSERT(q_min == q_max,Min/max values must match when n_q == 1)

     q(1) = q_min

     lam(1) = lambda(q(1), m, k)

  else

     $ABORT(n_q must be 1 or greater)

  endif

  ! Write out results

  hg = hgroup_t(filename, CREATE_FILE)
  call write_attr(hg, 'l', l)
  call write_attr(hg, 'm', m)
  call write_attr(hg, 'k', k)
  call write_dset(hg, 'rossby', rossby)
  call write_dset(hg, 'q', q)
  call write_dset(hg, 'lambda', lam)
  call hg%final()

  ! Finish

end program eval_lambda
