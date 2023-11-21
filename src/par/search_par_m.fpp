! Module  : search_par_m
! Purpose : mode search parameters
!
! Copyright 2020 Rich Townsend
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

module search_par_m

  ! Uses

  use kinds_m

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: search_par_t
     character(64) :: search_type = 'GRID'
  end type search_par_t

  ! Access specifiers

  private

  public :: search_par_t
  public :: read_search_par

  ! Procedures

contains

  subroutine read_search_par (unit, se_p)

    integer, intent(in)             :: unit
    type(search_par_t), intent(out) :: se_p

    character(LEN(se_p%search_type)) :: search_type

    namelist /search/ search_type

    ! Count the number of search namelists

    rewind(unit)

    n_se_p = 0

    count_loop : do
       read(unit, NML=search, END=100)
       n_se_p = n_se_p + 1
    end do count_loop

100 continue

    ! Read search parameters

    rewind(unit)

    ! Set default values

    se = search_par_t()

    search_type = se_p%search_type

    ! Read the namelist

    read(unit, NML=search)

    ! Store read values

    se_p%search_type = search_type

    ! Finish

    return

  end subroutine read_search_par

end module search_par_m
