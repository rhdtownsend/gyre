! Module   : gyre_out_par
! Purpose  : output parameters
!
! Copyright 2013-2018 Rich Townsend
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

module gyre_out_par

  ! Uses

  use core_kinds

  use gyre_constants

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: out_par_t
     character(256)          :: freq_units = 'NONE'
     character(256)          :: freq_frame = 'INERTIAL'
     character(FILENAME_LEN) :: summary_file = ''
     character(256)          :: summary_file_format = 'HDF'
     character(2048)         :: summary_item_list = 'l,n_pg,omega,freq'
     character(FILENAME_LEN) :: details_template = ''
     character(256)          :: details_file_format = 'HDF'
     character(2048)         :: details_item_list = 'l,n_pg,omega,freq,x,xi_r,xi_h'
     character(2048)         :: details_filter_list = ''
     character(256)          :: label = ''
     logical                 :: prune_details = .FALSE.
  end type out_par_t

 ! Access specifiers

  private

  public :: out_par_t
  public :: read_out_par

  ! Procedures

contains

  subroutine read_out_par (unit, stage, ot_p)

    integer, intent(in)          :: unit
    character(*), intent(in)     :: stage
    type(out_par_t), intent(out) :: ot_p

    integer                                  :: n_ot_p
    character(LEN(ot_p%freq_units))          :: freq_units
    character(LEN(ot_p%freq_frame))          :: freq_frame
    character(LEN(ot_p%summary_file))        :: summary_file
    character(LEN(ot_p%summary_file_format)) :: summary_file_format
    character(LEN(ot_p%summary_item_list))   :: summary_item_list
    character(LEN(ot_p%details_template))    :: details_template
    character(LEN(ot_p%details_file_format)) :: details_file_format
    character(LEN(ot_p%details_item_list))   :: details_item_list
    character(LEN(ot_p%details_filter_list)) :: details_filter_list
    character(LEN(ot_p%label))               :: label
    logical                                  :: prune_details

    namelist /ad_output/ freq_units, freq_frame, &
                         summary_file, summary_file_format, summary_item_list, &
                         details_template, details_file_format, details_item_list, &
                         label, prune_details

    namelist /nad_output/ freq_units, freq_frame, &
                          summary_file, summary_file_format, summary_item_list, &
                          details_template, details_file_format, details_item_list, details_filter_list, &
                          label, prune_details

    ! Count the number of output namelists

    rewind(unit)

    n_ot_p = 0

    count_loop : do
       select case (stage)
       case ('ad')
          read(unit, NML=ad_output, END=100)
       case ('nad')
          read(unit, NML=nad_output, END=100)
       case default
          $ABORT(Invalid stage)
       end select
       n_ot_p = n_ot_p + 1
    end do count_loop

100 continue

    $ASSERT(n_ot_p == 1,Input file should contain exactly one &ad_output and one &nad_output namelist)

    ! Read output parameters

    rewind(unit)

    ! Set default values

    ot_p = out_par_t()

    freq_units = ot_p%freq_units
    freq_frame = ot_p%freq_frame
    summary_file = ot_p%summary_file
    summary_file_format = ot_p%summary_file_format
    summary_item_list = ot_p%summary_item_list
    details_template = ot_p%details_template
    details_file_format = ot_p%details_file_format
    details_item_list = ot_p%details_item_list
    details_filter_list = ot_p%details_filter_list
    label = ot_p%label
    prune_details = ot_p%prune_details

    ! Read the namelist

    select case (stage)
    case ('ad')
       read(unit, NML=ad_output)
    case ('nad')
       read(unit, NML=nad_output)
    case default
       $ABORT(Invalid stage)
    end select

    ! Store read values

    ot_p%freq_units = freq_units
    ot_p%freq_frame = freq_frame
    ot_p%summary_file = summary_file
    ot_p%summary_file_format = summary_file_format
    ot_p%summary_item_list = summary_item_list
    ot_p%details_template = details_template
    ot_p%details_file_format = details_file_format
    ot_p%details_item_list = details_item_list
    ot_p%details_filter_list = details_filter_list
    ot_p%label = label
    ot_p%prune_details = prune_details

    ! Finish

    return

  end subroutine read_out_par

end module gyre_out_par
