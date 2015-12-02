! Module   : gyre_osc_par
! Purpose  : oscillation parameters
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

module gyre_osc_par

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: osc_par_t
     real(WP)        :: x_ref
     character(64)   :: rotation_method
     character(64)   :: variables_set
     character(64)   :: inner_bound
     character(64)   :: outer_bound
     character(64)   :: inertia_norm
     character(2048) :: tag_list
     logical         :: nonadiabatic
     logical         :: cowling_approx
     logical         :: narf_approx
     logical         :: reduce_order
  end type osc_par_t

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

  public :: osc_par_t
  public :: read_osc_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

!****

  subroutine read_osc_par (unit, os_p)

    integer, intent(in)                       :: unit
    type(osc_par_t), allocatable, intent(out) :: os_p(:)

    integer                              :: n_os_p
    integer                              :: i
    character(LEN(os_p%rotation_method)) :: rotation_method
    character(LEN(os_p%variables_set))   :: variables_set
    character(LEN(os_p%inner_bound))     :: inner_bound
    character(LEN(os_p%outer_bound))     :: outer_bound
    character(LEN(os_p%inertia_norm))    :: inertia_norm
    character(LEN(os_p%tag_list))        :: tag_list
    logical                              :: nonadiabatic
    logical                              :: cowling_approx
    logical                              :: narf_approx
    logical                              :: reduce_order
    real(WP)                             :: x_ref

    namelist /osc/ x_ref, rotation_method, inner_bound, outer_bound, variables_set, &
         inertia_norm, tag_list, nonadiabatic, cowling_approx, narf_approx, &
         reduce_order

    ! Count the number of osc namelists

    rewind(unit)

    n_os_p = 0

    count_loop : do
       read(unit, NML=osc, END=100)
       n_os_p = n_os_p + 1
    end do count_loop

100 continue

    ! Read oscillation parameters

    rewind(unit)

    allocate(os_p(n_os_p))

    read_loop : do i = 1,n_os_p

       x_ref = HUGE(0._WP)

       rotation_method = 'NULL'
       variables_set = 'DZIEM'
       inner_bound = 'REGULAR'
       outer_bound = 'ZERO'
       inertia_norm = 'BOTH'
       tag_list = ''

       nonadiabatic = .FALSE.
       cowling_approx = .FALSE.
       narf_approx = .FALSE.
       reduce_order = .TRUE.

       read(unit, NML=osc)

       ! Initialize the osc_par

       os_p(i) = osc_par_t(x_ref=x_ref, &
                           rotation_method=rotation_method, &
                           variables_set=variables_set, &
                           inner_bound=inner_bound, &
                           outer_bound=outer_bound, &
                           inertia_norm=inertia_norm, &
                           tag_list=tag_list, &
                           nonadiabatic=nonadiabatic, &
                           cowling_approx=cowling_approx, &
                           narf_approx=narf_approx, &
                           reduce_order=reduce_order)

    end do read_loop

    ! Finish

    return

  end subroutine read_osc_par

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (os_p, root_rank)

    type(osc_par_t), intent(inout) :: os_p$ARRAY_SPEC($RANK)
    integer, intent(in)            :: root_rank

    ! Broadcast the osc_par_t

    call bcast(os_p%x_ref, root_rank)

    call bcast(os_p%rotation_method, root_rank)
    call bcast(os_p%variables_set, root_rank)
    call bcast(os_p%inner_bound, root_rank)
    call bcast(os_p%outer_bound, root_rank)
    call bcast(os_p%inertia_norm, root_rank)
    call bcast(os_p%tag_list, root_rank)

    call bcast(os_p%nonadiabatic, root_rank)
    call bcast(os_p%cowling_approx, root_rank)
    call bcast(os_p%narf_approx, root_rank)
    call bcast(os_p%reduce_order, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(osc_par_t),0)
  $BCAST_ALLOC(type(osc_par_t),1)

  $endif

end module gyre_osc_par
