! Program  : gyre_trad
! Purpose  : traditional approximation support
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

module gyre_trad

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup
  use core_system

  use gyre_constants
  use gyre_trad_table

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Module variables

  logical, save            :: inited_m = .FALSE.
  type(trad_table_t), save :: tt_m

  ! Interfaces

  interface trad_lambda
     module procedure trad_lambda_r_
     module procedure trad_lambda_c_
  end interface trad_lambda

  ! Access specifiers

  private

  public :: trad_lambda

  ! Procedures

contains

  subroutine init_ ()

    type(hgroup_t) :: hg

    ! Load the trad_table_t from the data directory

    hg = hgroup_t(TRIM(GYRE_DIR)//'/data/trad_table.h5', OPEN_FILE)
    call read(hg, tt_m)
    call hg%final()

    ! Finish

    return

  end subroutine init_

  !****

  function trad_lambda_r_ (nu, l, m, rossby) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    logical, intent(in)  :: rossby
    real(WP)             :: lambda

    integer :: k
    
    ! Evaluate the eigenvalue of Laplace's tidal equation (real)

    if (.NOT. inited_m) then
       call init_()
       inited_m = .TRUE.
    endif

    $ASSERT(ABS(m) <= tt_m%m_max,Out-of-bounds of trad_table_t)

    if (rossby) then
       k = -(l - ABS(m) + 1)
    else
       k = l - ABS(m)
    endif

    $ASSERT(k >= tt_m%k_min,Out-of-bounds of trad_table_t)
    $ASSERT(k <= tt_m%k_max,Out-of-bounds of trad_table_t)

    $ASSERT(.NOT. (m == 0 .AND. k < 0),Invalid Rossby mode)

    if (m >= 0) then
       lambda = tt_m%tf(m,k)%lambda(nu)
    else
       lambda = tt_m%tf(-m,k)%lambda(-nu)
    endif

    ! Finish

    return

  end function trad_lambda_r_

  !****

  function trad_lambda_c_ (nu, l, m, rossby) result (lambda)

    complex(WP), intent(in) :: nu
    integer, intent(in)     :: l
    integer, intent(in)     :: m
    logical, intent(in)     :: rossby
    complex(WP)             :: lambda

    integer :: k

    ! Evaluate the eigenvalue of Laplace's tidal equation (complex)

    if (.NOT. inited_m) then
       call init_()
       inited_m = .TRUE.
    endif

    $ASSERT(ABS(m) <= tt_m%m_max,Out-of-bounds of trad_table_t)

    if (rossby) then
       k = -(l - ABS(m) + 1)
    else
       k = l - ABS(m)
    endif

    $ASSERT(k >= tt_m%k_min,Out-of-bounds of trad_table_t)
    $ASSERT(k <= tt_m%k_max,Out-of-bounds of trad_table_t)

    $ASSERT(.NOT. (m == 0 .AND. k < 0),Invalid Rossby mode)

    if (m >= 0) then
       lambda = tt_m%tf(m,k)%lambda(nu)
    else
       lambda = tt_m%tf(-m,k)%lambda(-nu)
    endif

    ! Finish

    return

  end function trad_lambda_c_

end module gyre_trad
