! Module   : gyre_twopt_file
! Purpose  : read two-point model files
!
! Copyright 2018 Rich Townsend
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

module gyre_twopt_file

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_model_par
  use gyre_twopt_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_twopt_model

  ! Procedures

contains

  subroutine read_twopt_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                      :: unit
    real(WP)                     :: V_o
    real(WP)                     :: dV_o
    real(WP)                     :: U_o
    real(WP)                     :: dU_o
    real(WP)                     :: c_1_i
    type(twopt_model_t), pointer :: tm

    namelist /model/ V_o, dV_o, U_o, dU_o, c_1_i

    ! Read data from the TWOPT-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from TWOPT file', TRIM(ml_p%file)
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif

    ! Read the file

    open(NEWUNIT=unit, file=ml_p%file, STATUS='OLD')

    read(unit, NML=model)

    close(unit)

    ! Initialize the model

    allocate(tm)

    tm = twopt_model_t(ml_p, V_o, dV_o, U_o, dU_o, c_1_i)

    ! Return a pointer

    ml => tm

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_twopt_model

end module gyre_twopt_file
