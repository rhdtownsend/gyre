! Module   : gyre_poly_file
! Purpose  : read POLY files
!
! Copyright 2013 Rich Townsend
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

module gyre_poly_file

  ! Uses

  use core_kinds
  use core_hgroup

  use gyre_mech_coeffs
  use gyre_poly_mech_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_poly_file

  ! Procedures

contains

  subroutine read_poly_file (file, deriv_type, mc, x)

    character(LEN=*), intent(in)                   :: file
    character(LEN=*), intent(in)                   :: deriv_type
    class(mech_coeffs_t), allocatable, intent(out) :: mc
    real(WP), allocatable, intent(out), optional   :: x(:)

    type(hgroup_t)        :: hg
    real(WP)              :: n_poly
    real(WP)              :: Gamma_1
    real(WP), allocatable :: xi(:)
    real(WP), allocatable :: Theta(:)
    real(WP), allocatable :: dTheta(:)

    ! Read the Lane-Emden solution from the POLY-format file

    call hg%init(file, OPEN_FILE)

    call read_attr(hg, 'n_poly', n_poly)
    call read_attr(hg, 'Gamma_1', Gamma_1)

    call read_dset_alloc(hg, 'xi', xi)
    call read_dset_alloc(hg, 'Theta', Theta)
    call read_dset_alloc(hg, 'dTheta', dTheta)

    call hg%final()

    ! Initialize the mech_coeffs

    allocate(poly_mech_coeffs_t::mc)

    select type (mc)
    type is (poly_mech_coeffs_t)
       call mc%init(xi,Theta, dTheta, n_poly, Gamma_1, deriv_type)
    class default
       $ABORT(Invalid mc type)
    end select

    ! If necessary, return the grid

    if(PRESENT(x)) then
       x = xi/xi(SIZE(xi))
    endif

    ! Finish

    return

  end subroutine read_poly_file

end module gyre_poly_file
