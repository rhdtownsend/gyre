! Program  : gyre_tar
! Purpose  : traditional approximation of rotation (TAR) eigenvalue evaluation
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_tar

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup
  use core_system

  use gyre_constants
  use gyre_tar_fit

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Module variables

  type(tar_fit_t), save :: tf_m
  logical, save         :: loaded_m = .FALSE.

  ! Interfaces

  interface tar_lambda
     module procedure tar_lambda_r_
     module procedure tar_lambda_c_
  end interface tar_lambda

  ! Access specifiers

  private

  public :: tar_lambda

  ! Procedures

contains

  subroutine load_fit_ (l, m, rossby)

    integer, intent(in) :: l
    integer, intent(in) :: m
    logical, intent(in) :: rossby

    integer                 :: k
    character(FILENAME_LEN) :: filename
    type(hgroup_t)          :: hg

    ! Load the appropriate tar_fit_t from the data directory

    if (rossby) then
       k = -(l - ABS(m) + 1)
    else
       k = l - ABS(m)
    endif

    if (loaded_m) then
       if (tf_m%m == m .AND. tf_m%k == k) return
    endif

    write(filename, 100) m, k
100 format(SP,'tar_fit.m',I0,'.k',I0,'.h5')

    hg = hgroup_t(TRIM(GYRE_DIR)//'/data/tar/'//TRIM(filename), OPEN_FILE)
    call read(hg, tf_m)
    call hg%final()

    loaded_m = .TRUE.

    ! Finish

    return

  end subroutine load_fit_

  !****

  $define $TAR_LAMBDA $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function tar_lambda_${SUFFIX}_ (nu, l, m, rossby) result (lambda)

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

  end function tar_lambda_${SUFFIX}_

  $endsub

  $TAR_LAMBDA(r,real)
  $TAR_LAMBDA(c,complex)

end module gyre_tar
