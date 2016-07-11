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
  use gyre_trad_fit

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Module variables

  type(trad_fit_t), save :: tf_m
  logical, save          :: loaded_m = .FALSE.

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

  subroutine load_fit_ (l, m, rossby)

    integer, intent(in) :: l
    integer, intent(in) :: m
    logical, intent(in) :: rossby

    integer                 :: k
    character(FILENAME_LEN) :: filename
    type(hgroup_t)          :: hg

    ! Load the appropriate trad_fit_t from the data directory

    if (rossby) then
       k = -(l - ABS(m) + 1)
    else
       k = l - ABS(m)
    endif

    if (loaded_m) then
       if (tf_m%m == m .AND. tf_m%k == k) return
    endif

    write(filename, 100) m, k
100 format(SP,'trad_fit.m',I0,'.k',I0,'.h5')

    hg = hgroup_t(TRIM(GYRE_DIR)//'/data/trad/'//TRIM(filename), OPEN_FILE)
    call read(hg, tf_m)
    call hg%final()

    loaded_m = .TRUE.

    ! Finish

    return

  end subroutine load_fit_

  !****

  $define $TRAD_LAMBDA $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function trad_lambda_${SUFFIX}_ (nu, l, m, rossby) result (lambda)

    $TYPE(WP), intent(in) :: nu
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    logical, intent(in)  :: rossby
    $TYPE(WP)            :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation

    call load_fit_(l, m, rossby)

    lambda = tf_m%lambda(nu)

    ! Finish

    return

  end function trad_lambda_${SUFFIX}_

  $endsub

  $TRAD_LAMBDA(r,real)
  $TRAD_LAMBDA(c,complex)

end module gyre_trad
