! Module   : gyre_oni_file
! Purpose  : read ONI files
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

module gyre_oni_file

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_model_par
  use gyre_oni_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_oni_model

  ! Procedures

contains

  subroutine read_oni_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                    :: unit
    character(LEN=2048)        :: line
    integer                    :: n
    real(WP), allocatable      :: M_r(:)
    real(WP), allocatable      :: r(:)
    integer                    :: j
    type(oni_model_t), pointer :: om

    ! Read the ONI-format file (ECSV)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from ONI file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    ! Open the file

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Count lines

    n = 0

    count_loop: do
       read(unit, 110, END=200) line
110     format(A)
       if (SCAN(line(1:1), '0123456789+-') /= 0) n = n+1
    end do count_loop
    
200 continue

    rewind(unit)

    allocate(M_r(n))
    allocate(r(n))

    read_loop: do j = 1, n
       do
          read(unit, 110) line
          if (SCAN(line(1:1), '0123456789+-') /= 0) then
             read(line, *) M_r(j), r(j)
             exit
          end if
       end do
    end do read_loop
    
    close(unit)

    ! Initialize the oni_model_t

    allocate(om, SOURCE=oni_model_t(M_r, r, ml_p%Gamma_1))

    ! Return a pointer

    ml => om

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_oni_model

end module gyre_oni_file
