! Module  : parfait_file_m
! Purpose : read PARFAIT model files
!
! Copyright 2022-2023 Rich Townsend & The GYRE Team
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

module parfait_file_m

  ! Uses

  use kinds_m
  use hgroup_m

  use model_m
  use model_par_m
  use parfait_model_m
  use util_m

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

    real(WP)                       :: y_c
    real(WP)                       :: z_s
    real(WP), allocatable          :: x(:)
    real(WP), allocatable          :: d(:)
    real(WP), allocatable          :: Gamma_1(:)
    type(parfait_model_t), pointer :: pm

    ! Read data from the PARFAIT-format file

    call read_parfait_data(ml_p%file, y_c, z_s, x, d, Gamma_1)

    ! Initialize the parfait_model_t

    allocate(pm, SOURCE=parfait_model_t(x, d, Gamma_1, y_c, z_s, ml_p))

    ! Return a pointer

    ml => pm

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_parfait_model

  !****

  subroutine read_parfait_data (file, y_c, z_s, x, d, Gamma_1)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: y_c
    real(WP), intent(out)              :: z_s
    real(WP), allocatable, intent(out) :: x(:)
    real(WP), allocatable, intent(out) :: d(:)
    real(WP), allocatable, intent(out) :: Gamma_1(:)

    type(hgroup_t) :: hg
    character(16)  :: type
    integer        :: version
    integer        :: N

    ! Read data from the PARFAIT-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from PARFAIT file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(file)
110    format(3X,A,1X,A)
    endif

    hg = hgroup_t(file, OPEN_FILE_RO)

    ! Read the header

    call read_attr(hg, 'type', type)

    if (type /= 'PARFAIT') then
       $ABORT(File type mismatch)
    end if
    
    call read_attr(hg, 'version', version)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'File version', version/100._WP
120    format(3X,A,1X,F4.2,1X,A)
    endif

    call read_attr(hg, 'N', N)

    ! Read the data

    select case (version)
    case (100)
       call read_parfait_data_v1_00_()
    case default
       $ABORT(Unrecognized PARFAIT file version)
    end select

    call hg%final()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', N, 'shells'
130    format(3X,A,1X,I0,1X,A)
    endif

    ! Finish

    return

  contains

    subroutine read_parfait_data_v1_00_ ()

      ! Read data from the version-1.00 file

      call read_attr(hg, 'y_c', y_c)
      call read_attr(hg, 'z_s', z_s)

      allocate(x(N+1))
      allocate(d(N))
      allocate(Gamma_1(N))

      call read_dset(hg, 'x', x)
      call read_dset(hg, 'd', d)
      call read_dset(hg, 'Gamma_1', Gamma_1)

      ! Finish

      return

    end subroutine read_parfait_data_v1_00_

  end subroutine read_parfait_data

end module parfait_file_m
