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

  type(tar_fit_t), save, volatile :: tf_m
  logical, save, volatile         :: loaded_m = .FALSE.

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

  $define $TAR_LAMBDA $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function tar_lambda_${SUFFIX}_ (nu, l, m, rossby) result (lambda)

    $TYPE(WP), intent(in) :: nu
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    logical, intent(in)  :: rossby
    $TYPE(WP)            :: lambda

    integer                 :: k
    character(FILENAME_LEN) :: filename
    type(hgroup_t)          :: hg
    
    ! Evaluate the eigenvalue of Laplace's tidal equation

    ! Set up k

    if (rossby) then
       k = -(l - ABS(m) + 1)
    else
       k = l - ABS(m)
    endif

    ! If necessary, load the appropriate tf_m from the data directory.
    ! Use a double-checked locking pattern (DCLP) to avoid the overhead
    ! of the CRITICAL section (although see 'C++ and the Perils of
    ! Double-Checked Locking', Meyers & Alexandrescu 2004; I think the
    ! use of volatile in the declarations of tf_m and loaded_m should
    ! make things work)

    if (must_load_(k, m)) then

       !$OMP CRITICAL

       if (must_load_(k, m)) then

          write(filename, 100) m, k
100       format(SP,'tar_fit.m',I0,'.k',I0,'.h5')

          hg = hgroup_t(TRIM(GYRE_DIR)//'/data/tar/'//TRIM(filename), OPEN_FILE)
          call read(hg, tf_m)
          call hg%final()
       
          loaded_m = .TRUE.

       endif

       !$OMP END CRITICAL

    endif

    ! The following assertion will fail if this routine is called
    ! inside a parallel section with each thread having a different l,
    ! m, or rossby

    $ASSERT_DEBUG(tf_m%k == k .AND. tf_m%m == m,tar_fit has wrong parameters)

    lambda = tf_m%lambda(nu)

    ! Finish

    return

  contains

    function must_load_ (k, m) result (must_load)

      integer, intent(in) :: k
      integer, intent(in) :: m
      logical             :: must_load

      ! Decide whether tf_m must be loaded

      if (loaded_m) then
         must_load = tf_m%k /= k .OR. tf_m%m /= m
      else
         must_load = .TRUE.
      end if

      ! Finish

      return

    end function must_load_

  end function tar_lambda_${SUFFIX}_

  $endsub

  $TAR_LAMBDA(r,real)
  $TAR_LAMBDA(c,complex)

end module gyre_tar
