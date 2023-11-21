! Module  : poly_file_m
! Purpose : read POLY files
!
! Copyright 2013-2022 Rich Townsend & The MESA Team
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

module poly_file_m

  ! Uses

  use kinds_m
  use hgroup_m

  use model_m
  use model_par_m
  use poly_model_m
  use util_m

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
    integer                     :: n_r
    real(WP), allocatable       :: Delta_b(:)
    real(WP)                    :: Gamma_1
    real(WP), allocatable       :: z(:)
    real(WP), allocatable       :: theta(:)
    real(WP), allocatable       :: dtheta(:)
    real(WP)                    :: Omega_rot
    type(poly_model_t), pointer :: pm

    ! Read the POLY-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from POLY file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    hg = hgroup_t(ml_p%file, OPEN_FILE_RO)

    call read_attr(hg, 'n_r', n_r)

    call read_attr_alloc(hg, 'n_poly', n_poly)
    if (n_r > 1) then
       call read_attr_alloc(hg, 'Delta_b', Delta_b)
    else
       allocate(Delta_b(0))
    endif
    call read_attr(hg, 'Gamma_1', Gamma_1)

    call read_dset_alloc(hg, 'z', z)
    call read_dset_alloc(hg, 'theta', theta)
    call read_dset_alloc(hg, 'dtheta', dtheta)

    call hg%final()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Read', SIZE(z), 'points'
110    format(3X,A,1X,I0,1X,A)
    endif

    ! Set up Omega_rot

    Omega_rot = 0._WP

    ! Initialize the poly_model_t

    allocate(pm, SOURCE=poly_model_t(z, theta, dtheta, n_poly, Delta_b, Gamma_1, Omega_rot))

    ! Return a pointer

    ml => pm

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_poly_model

end module poly_file_m
