! Program  : gyre_trad
! Purpose  : traditional approximation support
!
! Copyright 2013-2014 Rich Townsend
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

  use gyre_trad_table

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Module variables

  type(trad_table_t), save :: tt_m

  ! Interfaces

  interface trad_lambda
     module procedure trad_lambda_r_
     module procedure trad_lambda_c_
  end interface trad_lambda

  ! Access specifiers

  private

  public :: init_trad
  public :: trad_lambda

  ! Procedures

contains

  subroutine init_trad (gyre_dir)

    character(*), intent(in) :: gyre_dir

    type(hgroup_t) :: hg

    ! Load the trad_table_t from the data directory

    hg = hgroup_t(TRIM(gyre_dir)//'/data/trad_table.h5', OPEN_FILE)
    call read(hg, tt_m)
    call hg%final()

    ! Finish

    return

  end subroutine init_trad

!****

  function trad_lambda_r_ (nu, l, m) result (lambda)

    integer, intent(in)  :: l
    integer, intent(in)  :: m
    real(WP), intent(in) :: nu
    real(WP)             :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (real)

    lambda = tt_m%tf(l,m)%lambda(nu)

    ! Finish

    return

  end function trad_lambda_r_

!****

  function trad_lambda_c_ (nu, l, m) result (lambda)

    integer, intent(in)     :: l
    integer, intent(in)     :: m
    complex(WP), intent(in) :: nu
    complex(WP)             :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (complex)

    lambda = tt_m%tf(l,m)%lambda(nu)

    ! Finish

    return

  end function trad_lambda_c_

end module gyre_trad
