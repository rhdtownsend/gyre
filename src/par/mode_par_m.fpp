! Module  : mode_par_m
! Purpose : mode parameters
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module mode_par_m

  ! Uses

  use kinds_m

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_par_t
     integer       :: l = 0
     integer       :: m = 0
     integer       :: n_pg_min = -HUGE(0)
     integer       :: n_pg_max = HUGE(0)
     character(64) :: tag = ''
  end type mode_par_t

 ! Access specifiers

  private

  public :: mode_par_t
  public :: read_mode_par

  ! Procedures

contains

  subroutine read_mode_par (unit, md_p)

    integer, intent(in)                        :: unit
    type(mode_par_t), allocatable, intent(out) :: md_p(:)

    integer                  :: n_md_p
    integer                  :: i
    integer                  :: l
    integer                  :: m
    integer                  :: n_pg_min
    integer                  :: n_pg_max
    character(LEN(md_p%tag)) :: tag
 
    namelist /mode/ l, m, n_pg_min, n_pg_max, tag

    ! Count the number of mode namelists

    rewind(unit)

    n_md_p = 0

    count_loop : do
       read(unit, NML=mode, END=100)
       n_md_p = n_md_p + 1
    end do count_loop

100 continue

    ! Read mode parameters

    rewind(unit)

    allocate(md_p(n_md_p))

    read_loop : do i = 1, n_md_p

       ! Set default values

       md_p(i) = mode_par_t()

       l = md_p(i)%l
       m = md_p(i)%m
       n_pg_min = md_p(i)%n_pg_min
       n_pg_max = md_p(i)%n_pg_max
       tag = md_p(i)%tag

       ! Read the namelist

       read(unit, NML=mode)

       ! Store read values

       md_p(i)%l = l
       md_p(i)%m = m
       md_p(i)%n_pg_min = n_pg_min
       md_p(i)%n_pg_max = n_pg_max
       md_p(i)%tag = tag

    end do read_loop

    ! Finish

    return

  end subroutine read_mode_par

end module mode_par_m
