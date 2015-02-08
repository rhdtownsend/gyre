! Module   : gyre_out_par
! Purpose  : output parameters
!
! Copyright 2013-2015 Rich Townsend
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
$include 'core_parallel.inc'

module gyre_out_par

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: out_par_t
     character(256)          :: freq_units
     character(256)          :: freq_frame
     character(FILENAME_LEN) :: summary_file
     character(256)          :: summary_file_format
     character(2048)         :: summary_item_list
     character(FILENAME_LEN) :: mode_prefix
     character(FILENAME_LEN) :: mode_template
     character(256)          :: mode_file_format
     character(2048)         :: mode_item_list
     character(256)          :: label
     logical                 :: prune_modes
  end type out_par_t

  ! Interfaces

  $if ($MPI)

  interface bcast
     module procedure bcast_0_
     module procedure bcast_1_
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_0_
     module procedure bcast_alloc_1_
  end interface bcast_alloc

  $endif

 ! Access specifiers

  private

  public :: out_par_t
  public :: read_out_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  subroutine read_out_par (unit, up)

    integer, intent(in)          :: unit
    type(out_par_t), intent(out) :: up

    integer                                :: n_up
    character(LEN(up%freq_units))          :: freq_units
    character(LEN(up%freq_frame))          :: freq_frame
    character(LEN(up%summary_file))        :: summary_file
    character(LEN(up%summary_file_format)) :: summary_file_format
    character(LEN(up%summary_item_list))   :: summary_item_list
    character(LEN(up%mode_prefix))         :: mode_prefix
    character(LEN(up%mode_template))       :: mode_template
    character(LEN(up%mode_file_format))    :: mode_file_format
    character(LEN(up%mode_item_list))      :: mode_item_list
    character(LEN(up%label))               :: label
    logical                                :: prune_modes

    namelist /output/ freq_units, freq_frame, summary_file, summary_file_format, summary_item_list, &
                      mode_prefix, mode_template, mode_file_format, mode_item_list, prune_modes

    ! Count the number of output namelists

    rewind(unit)

    n_up = 0

    count_loop : do
       read(unit, NML=output, END=100)
       n_up = n_up + 1
    end do count_loop

100 continue

    $ASSERT(n_up == 1,Input file should contain exactly one &output namelist)

    ! Read output parameters

    freq_units = 'NONE'
    freq_frame = 'INERTIAL'

    summary_file = ''
    summary_file_format = 'HDF'
    summary_item_list = 'l,n_pg,omega,freq'
    
    mode_prefix = ''
    mode_template = ''
    mode_file_format = 'HDF'
    mode_item_list = TRIM(summary_item_list)//',x,xi_r,xi_h'

    label = ''

    prune_modes = .FALSE.

    rewind(unit)
    read(unit, NML=output)

    ! Initialize the out_par

    up = out_par_t(freq_units=freq_units, &
                  freq_frame=freq_frame, &
                  summary_file=summary_file, &
                  summary_file_format=summary_file_format, &
                  summary_item_list=summary_item_list, &
                  mode_prefix=mode_prefix, &
                  mode_template=mode_template, &
                  mode_file_format=mode_file_format, &
                  mode_item_list=mode_item_list, &
                  label=label, &
                  prune_modes=prune_modes)

    ! Finish

    return

  end subroutine read_out_par

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (up, root_rank)

    type(out_par_t), intent(inout) :: up$ARRAY_SPEC($RANK)
    integer, intent(in)            :: root_rank

    ! Broadcast the out_par_t

    call bcast(up%freq_units, root_rank)
    call bcast(up%freq_frame, root_rank)
    call bcast(up%summary_file, root_rank)
    call bcast(up%summary_file_format, root_rank)
    call bcast(up%summary_item_list, root_rank)
    call bcast(up%mode_prefix, root_rank)
    call bcast(up%mode_template, root_rank)
    call bcast(up%mode_file_format, root_rank)
    call bcast(up%mode_item_list, root_rank)
    call bcast(up%label, root_rank)
    
    call bcast(up%prune_modes, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(out_par_t),0)
  $BCAST_ALLOC(type(out_par_t),1)

  $endif

end module gyre_out_par
