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
  integer                   :: k
  real(WP)                  :: q_min
  real(WP)                  :: q_max
  integer                   :: n_q
  logical                   :: log_q
  character(:), allocatable :: filename

  real(WP), allocatable :: q(:)
  real(WP), allocatable :: lam(:)
  integer               :: i
  type(hgroup_t)        :: hg

  ! Read parameters

  $ASSERT(n_arg() == 7,Syntax: eval_lambda m k q_min q_max n_q log_q filename)

  call get_arg(1, m)
  call get_arg(2, k)
  call get_arg(3, q_min)
  call get_arg(4, q_max)
  call get_arg(5, n_q)
  call get_arg(6, log_q)
  call get_arg(7, filename)

  ! Initialize

  call init_parallel()

  call init_math()

  ! Allocate arrays

  allocate(q(n_q))
  allocate(lam(n_q))

  ! Evaluate them

  !$OMP PARALLEL DO SCHEDULE (GUIDED)
  do i = 1, n_q

     if (log_q) then
        q(i) = 10**((LOG10(q_min)*(n_q-i) + LOG10(q_max)*(i-1))/(n_q-1))
     else
        q(i) = (q_min*(n_q-i) + q_max*(i-1))/(n_q-1)
     endif

     lam(i) = lambda(q(i), m, k)

  end do

  ! Write out results

  hg = hgroup_t(filename, CREATE_FILE)
  call write_attr(hg, 'm', m)
  call write_attr(hg, 'k', k)
  call write_dset(hg, 'q', q)
  call write_dset(hg, 'lambda', lam)
  call hg%final()

  ! Finish

end program eval_lambda
