! Program  : test_trad_table
! Purpose  : test trad_table_t types
!
! Copyright 2013-2015 Rich Townsend
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

program test_trad_table

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system

  use astro_hough

  use gyre_trad_func
  use gyre_trad_table

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable :: filename
  real(WP)                  :: nu_min
  real(WP)                  :: nu_max
  integer                   :: n_nu
  real(WP)                  :: lambda_tol

  type(hgroup_t)        :: hg
  type(trad_table_t)    :: tt
  integer               :: m
  integer               :: k
  integer               :: i
  real(WP)              :: nu
  real(WP)              :: lambda_tab
  real(WP)              :: lambda_chk
  real(WP)              :: lambda_nrm
  real(WP), allocatable :: err(:)

  ! Read parameters

  $ASSERT(n_arg() == 5,Syntax: test_trad_table filename nu_min nu_max n_nu lambda_tol)

  call get_arg(1, filename)
  call get_arg(2, nu_min)
  call get_arg(3, nu_max)
  call get_arg(4, n_nu)
  call get_arg(5, lambda_tol)

  ! Read the trad_table_t

  hg = hgroup_t(filename, OPEN_FILE)
  call read(hg, tt)
  call hg%final()

  ! Test the entries

  allocate(err(n_nu))

  do m = 0, tt%m_max
     do k = tt%k_min, tt%k_max

        if (m == 0 .AND. k < 0) cycle

        do i = 1, n_nu

           nu = (nu_min*(n_nu-i) + nu_max*(i-1))/(n_nu-1)

           lambda_tab = tt%tf(m,k)%lambda(nu)
           lambda_chk = lambda(nu, m, k, lambda_tol)

           lambda_nrm = MIN(ABS(lambda_tab), ABS(lambda_chk))

           if (lambda_nrm /= 0._WP) then
              err(i) = (lambda_tab - lambda_chk)/lambda_nrm
           else
              err(i) = lambda_tab - lambda_chk
           endif

        end do

        write(OUTPUT_UNIT, *) 'Processed:', m, k, MINVAL(err), MAXVAL(err)

     end do
  end do

  ! Finish

end program test_trad_table
