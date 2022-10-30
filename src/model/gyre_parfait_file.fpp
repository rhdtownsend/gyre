! Module   : gyre_parfait_file
! Purpose  : read PARFAIT model files
!
! Copyright 2022 Rich Townsend & The GYRE Team
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

module gyre_parfait_file

  ! Uses

  use core_kinds
  use core_hgroup

  use gyre_model
  use gyre_model_par
  use gyre_parfait_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_parfait_model

  ! Procedures

contains

  subroutine read_parfait_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    type(hgroup_t)                 :: hg
    real(WP)                       :: y_c
    real(WP)                       :: z_s
    real(WP), allocatable          :: x(:)
    real(WP), allocatable          :: d(:)
    real(WP), allocatable          :: Gamma_1(:)
    type(parfait_model_t), pointer :: pm

    ! Read the PARFAIT-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from PARFAIT file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    hg = hgroup_t(ml_p%file, OPEN_FILE_RO)

    call read_attr(hg, 'y_c', y_c)
    call read_attr(hg, 'z_s', z_s)

    call read_dset_alloc(hg, 'x', x)
    call read_dset_alloc(hg, 'd', d)

    call read_dset_alloc(hg, 'Gamma_1', Gamma_1)

    call hg%final()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Read', SIZE(x), 'points'
110    format(3X,A,1X,I0,1X,A)
    endif

    ! Initialize the parfait_model_t

    allocate(pm, SOURCE=parfait_model_t(x, d, Gamma_1, y_c, z_s))

    ! Return a pointer

    ml => pm

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_parfait_model

end module gyre_parfait_file
