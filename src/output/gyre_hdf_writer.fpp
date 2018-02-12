! Module   : gyre_hdf_writer
! Purpose  : write HDF data
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

module gyre_hdf_writer

  ! Uses

  use core_kinds
  use core_hgroup

  use gyre_writer

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (writer_t) :: hdf_writer_t
     private
     type(hgroup_t) :: hg
   contains
     private
     procedure, public :: final
     procedure, public :: write_i_0
     procedure, public :: write_i_1
     procedure, public :: write_r_0
     procedure, public :: write_r_1
     procedure, public :: write_c_0
     procedure, public :: write_c_1
     procedure, public :: write_a_0
     procedure, public :: write_a_1
  end type hdf_writer_t

  ! Interfaces

  interface hdf_writer_t
     module procedure hdf_writer_t_
  end interface hdf_writer_t

  ! Access specifiers

  private

  public :: hdf_writer_t

  ! Procedures

contains

  function hdf_writer_t_ (file_name, label) result (wr)

    character(*), intent(in) :: file_name
    character(*), intent(in) :: label
    type(hdf_writer_t)       :: wr

    ! Construct the hdf_writer_t

    wr%hg = hgroup_t(file_name, CREATE_FILE)

    call write_attr(wr%hg, 'label', label)

    ! Finish

    return

  end function hdf_writer_t_

!****

  subroutine final (this)

    class(hdf_writer_t), intent(inout) :: this

    ! Finalize the hdf_writer

    call this%hg%final()

    ! Finish

    return

  end subroutine final

!****
  
  $define $WRITE $sub

  $local $INFIX $1
  $local $DATA_TYPE $2
  $local $DATA_RANK $3

  subroutine write_${INFIX}_${DATA_RANK} (this, name, data)

    class(hdf_writer_t), intent(inout) :: this
    character(*), intent(in)           :: name
    $DATA_TYPE, intent(in)             :: data$ARRAY_SPEC($DATA_RANK)

    ! Write the data

    $if($DATA_RANK == 1)
    call write_dset(this%hg, name, data)
    $else
    call write_attr(this%hg, name, data)
    $endif

    ! Finish

    return

  end subroutine write_${INFIX}_${DATA_RANK}

  $endsub

  $WRITE(i,integer,0)
  $WRITE(i,integer,1)

  $WRITE(r,real(WP),0)
  $WRITE(r,real(WP),1)

  $WRITE(c,complex(WP),0)
  $WRITE(c,complex(WP),1)

  $WRITE(a,character(*),0)
  $WRITE(a,character(*),1)

end module gyre_hdf_writer
