! Program  : gyre_bin
! Purpose  : binary orbital evolution code
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

program gyre_bin

  ! Uses

  use core_kinds, only : WP
  use core_parallel
  use core_system

  use gyre_constants
  use gyre_grid_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_tide
  use gyre_tide_par
  use gyre_util
  use gyre_version
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable      :: filename
  integer                        :: unit
  type(model_par_t)              :: ml_p
  type(osc_par_t), allocatable   :: os_p(:)
  type(num_par_t), allocatable   :: nm_p(:)
  type(grid_par_t), allocatable  :: gr_p(:)
  type(tide_par_t), allocatable  :: td_p(:)
  class(model_t), pointer        :: ml => null()
  integer                        :: i
  real(WP)                       :: tau_tot

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_bin <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_bin ['//VERSION//']', '-')
100  format(A)

     if (check_log_level('DEBUG')) then
        write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
        write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110     format(A,1X,A)
     endif

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

     write(OUTPUT_UNIT, *)

  endif

  ! Read the namelist file

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_tide_par(unit, td_p)

  close(unit)

  $ASSERT(SIZE(os_p) == 1,Must be exactly one mode parameter)
  $ASSERT(SIZE(nm_p) == 1,Must be exactly one num parameter)
  $ASSERT(SIZE(gr_p) == 1,Must be exactly one grid parameter)
  $ASSERT(SIZE(td_p) >= 1,Must be at least one tide parameter)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '-')
  endif

  ml => model_t(ml_p)

  ! Loop over tide parameters

  tide_par_loop : do i = 1, SIZE(td_p)

     ! Evaluate the tide

     tau_tot = 0._WP

     call eval_tide(ml, process_wave, os_p(1), nm_p(1), gr_p(1), td_p(i))

     print *,tau_tot

  end do tide_par_loop

  ! Clean up

  deallocate(ml)

  ! Finish

  call final_parallel()

contains

  subroutine process_wave (wv)

    type(wave_t), intent(in) :: wv

    ! Accumulate the torque

    print *,'Total torque:', wv%l, wv%m, wv%omega, wv%tau_ss()

!    tau_tot = tau_tot + wv%tau_ss()

    ! Finish

    return

  end subroutine process_wave

end program gyre_bin
