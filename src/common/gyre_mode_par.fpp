! Module   : gyre_mode_par
! Purpose  : mode parameters
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

module gyre_mode_par

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_par_t
     integer       :: l
     integer       :: m
     integer       :: n_pg_min
     integer       :: n_pg_max
     character(64) :: tag
  end type mode_par_t

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

  public :: mode_par_t
  public :: read_mode_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

!****

  subroutine read_mode_par (unit, mp)

    integer, intent(in)                        :: unit
    type(mode_par_t), allocatable, intent(out) :: mp(:)

    integer                :: n_mp
    integer                :: i
    integer                :: l
    integer                :: m
    integer                :: n_pg_min
    integer                :: n_pg_max
    character(LEN(mp%tag)) :: tag

    namelist /mode/ l, m, n_pg_min, n_pg_max, tag

    ! Count the number of mode namelists

    rewind(unit)

    n_mp = 0

    count_loop : do
       read(unit, NML=mode, END=100)
       n_mp = n_mp + 1
    end do count_loop

100 continue

    ! Read mode parameters

    rewind(unit)

    allocate(mp(n_mp))

    read_loop : do i = 1,n_mp

       l = 0
       m = 0

       n_pg_min = -HUGE(0)
       n_pg_max = HUGE(0)

       tag = ''

       read(unit, NML=mode)

       ! Initialize the mode_par

       mp(i) = mode_par_t(l=l, m=m, n_pg_min=n_pg_min, n_pg_max=n_pg_max, tag=tag)

    end do read_loop

    ! Finish

    return

  end subroutine read_mode_par

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (op, root_rank)

    type(mode_par_t), intent(inout) :: op$ARRAY_SPEC($RANK)
    integer, intent(in)             :: root_rank

    ! Broadcast the mode_par_t

    call bcast(op%l, root_rank)
    call bcast(op%m, root_rank)

    call bcast(op%n_pg_min, root_rank)
    call bcast(op%n_pg_max, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(mode_par_t),0)
  $BCAST_ALLOC(type(mode_par_t),1)

  $endif

end module gyre_mode_par
