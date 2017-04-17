! Module   : gyre_poly_file
! Purpose  : read POLY files
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

module gyre_poly_file

  ! Uses

  use core_kinds
  use core_hgroup

  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_poly_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_poly_model

  ! Procedures

contains

  subroutine read_poly_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    type(hgroup_t)              :: hg
    real(WP), allocatable       :: n_poly(:)
    integer                     :: n_d
    real(WP), allocatable       :: Delta_d(:)
    real(WP)                    :: Gamma_1
    real(WP), allocatable       :: xi(:)
    real(WP), allocatable       :: Theta(:)
    real(WP), allocatable       :: dTheta(:)
    real(WP)                    :: Omega_rot
    type(poly_model_t), pointer :: pm

    ! Read the POLY-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from POLY file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    hg = hgroup_t(ml_p%file, OPEN_FILE)

    call read_attr(hg, 'n_d', n_d)

    call read_attr_alloc(hg, 'n_poly', n_poly)
    if (n_d > 0) then
       call read_attr_alloc(hg, 'Delta_d', Delta_d)
    else
       allocate(Delta_d(0))
    endif
    call read_attr(hg, 'Gamma_1', Gamma_1)

    call read_dset_alloc(hg, 'xi', xi)
    call read_dset_alloc(hg, 'Theta', Theta)
    call read_dset_alloc(hg, 'dTheta', dTheta)

    call hg%final()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Read', SIZE(xi), 'points'
110    format(3X,A,1X,I0,1X,A)
    endif

    ! Set up Omega_rot

    if (ml_p%uniform_rot) then
       Omega_rot = uniform_Omega_rot(ml_p)
    else
       Omega_rot = 0._WP
    endif

    ! Initialize the poly_model_t

    allocate(pm, SOURCE=poly_model_t(xi, Theta, dTheta, n_poly, Delta_d, Gamma_1, Omega_rot))

    ! Return a pointer

    ml => pm

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_poly_model

end module gyre_poly_file
