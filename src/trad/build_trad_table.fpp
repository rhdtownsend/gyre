! Program  : build_trad_table
! Purpose  : build trad_table_t types
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

program build_trad_table

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system

  use gyre_trad_func
  use gyre_trad_table

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                   :: m_max
  integer                   :: k_min
  integer                   :: k_max
  real(WP)                  :: lambda_tol
  real(WP)                  :: cheby_tol
  character(:), allocatable :: filename

  type(trad_table_t) :: tt
  integer            :: m
  integer            :: k
  type(hgroup_t)     :: hg

  ! Read parameters

  $ASSERT(n_arg() == 6,Syntax: build_trad_table m_max k_min k_max lambda_tol cheby_tol filename)

  call get_arg(1, m_max)
  call get_arg(2, k_min)
  call get_arg(3, k_max)
  call get_arg(4, lambda_tol)
  call get_arg(5, cheby_tol)
  call get_arg(6, filename)

  ! Construct the trad_table_t

  tt = trad_table_t(m_max, k_min, k_max)

  ! Fill in the entries

  do m = 0, m_max
     do k = k_min, k_max
        write(OUTPUT_UNIT, *) 'Processing:', m, k
        if (m == 0 .AND. k < 0) cycle
        tt%tf(m,k) = trad_func_t(m, k, lambda_tol, cheby_tol)
     end do
  end do

  ! Write it out

  hg = hgroup_t(filename, CREATE_FILE)
  call write(hg, tt)
  call hg%final()

  ! Finish

end program build_trad_table
