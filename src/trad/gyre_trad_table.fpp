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
     integer                        :: m_max
     integer                        :: k_min
     integer                        :: k_max
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

  function trad_table_t_ (m_max, k_min, k_max) result (tt)

    integer, intent(in) :: m_max 
    integer, intent(in) :: k_min
    integer, intent(in) :: k_max 
    type(trad_table_t)  :: tt

    ! Construct the trad_table_t

    allocate(tt%tf(0:m_max,k_min:k_max))

    tt%m_max = m_max

    tt%k_min = k_min
    tt%k_max = k_max

    ! Finish

    return

  end function trad_table_t_

  !****

  $if ($HDF5)

  subroutine read_ (hg, tt)

    type(hgroup_t), intent(inout)   :: hg
    type(trad_table_t), intent(out) :: tt

    integer        :: m_max
    integer        :: k_min
    integer        :: k_max
    integer        :: m
    integer        :: k
    type(hgroup_t) :: hg_comp

    ! Read the trad_table_t

    call read_attr(hg, 'm_max', m_max)

    call read_attr(hg, 'k_min', k_min)
    call read_attr(hg, 'k_max', k_max)

    tt = trad_table_t(m_max, k_min, k_max)

    do m = 0, m_max
       do k = k_min, k_max
          if (m == 0 .AND. k < 0) cycle
          hg_comp = hgroup_t(hg, elem_group_name('tf', [m,k]))
          call read(hg_comp, tt%tf(m,k))
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

    integer        :: m
    integer        :: k
    type(hgroup_t) :: hg_comp

    ! Write the trad_table_t

    call write_attr(hg, 'm_max', tt%m_max)

    call write_attr(hg, 'k_min', tt%k_min)
    call write_attr(hg, 'k_max', tt%k_max)

    do m = 0, tt%m_max
       do k = tt%k_min, tt%k_max
          if (m == 0 .AND. k < 0) cycle
          hg_comp = hgroup_t(hg, elem_group_name('tf', [m,k]))
          call write(hg_comp, tt%tf(m,k))
          call hg_comp%final()
       end do
    end do

    ! Finish

    return

  end subroutine write_

  $endif

end module gyre_trad_table
