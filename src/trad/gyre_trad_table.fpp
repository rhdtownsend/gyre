! Program  : gyre_trad_table
! Purpose  : tables of trad_func_t types
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

module gyre_trad_table

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel

  use gyre_trad_func

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: trad_table_t
     type(trad_func_t), allocatable :: tf(:,:)
     integer                        :: l_max
  end type trad_table_t

  ! Interfaces

  interface trad_table_t
     module procedure trad_table_t_
  end interface trad_table_t

  $if ($HDF5)
  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write
  $endif

  ! Access specifiers

  private

  public :: trad_table_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif

  ! Procedures

contains

  function trad_table_t_ (l_max) result (tt)

    integer, intent(in) :: l_max 
    type(trad_table_t)  :: tt

    ! Construct the trad_table_t

    allocate(tt%tf(0:l_max,-l_max:l_max))

    tt%l_max = l_max

    ! Finish

    return

  end function trad_table_t_

!****

  $if ($HDF5)

  subroutine read_ (hg, tt)

    type(hgroup_t), intent(inout)   :: hg
    type(trad_table_t), intent(out) :: tt

    integer        :: l_max
    integer        :: l
    integer        :: m
    type(hgroup_t) :: hg_comp

    ! Read the trad_table_t

    call read_attr(hg, 'l_max', l_max)

    tt = trad_table_t(l_max)

    do l = 0, l_max
       do m = -l, l
          hg_comp = hgroup_t(hg, elem_group_name('tf', [l,m]))
          call read(hg_comp, tt%tf(l,m))
          call hg_comp%final()
       end do
    end do

    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, tt)

    type(hgroup_t), intent(inout)  :: hg
    type(trad_table_t), intent(in) :: tt

    integer        :: l
    integer        :: m
    type(hgroup_t) :: hg_comp

    ! Write the trad_table_t

    call write_attr(hg, 'l_max', tt%l_max)

    do l = 0, tt%l_max
       do m = -l, l
          hg_comp = hgroup_t(hg, elem_group_name('tf', [l,m]))
          call write(hg_comp, tt%tf(l,m))
          call hg_comp%final()
       end do
    end do

    ! Finish

    return

  end subroutine write_

  $endif

end module gyre_trad_table
