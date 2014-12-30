! Module   : gyre_oscpar
! Purpose  : oscillation parameters
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_oscpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: oscpar_t
     real(WP)        :: x_ref
     character(64)   :: rotation_method
     character(64)   :: variables_type
     character(64)   :: inner_bound_type
     character(64)   :: outer_bound_type
     character(64)   :: inertia_norm_type
     character(2048) :: tag_list
     logical         :: nonadiabatic
     logical         :: reduce_order
  end type oscpar_t

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

  public :: oscpar_t
  public :: read_oscpar
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

!****

  subroutine read_oscpar (unit, op)

    integer, intent(in)                      :: unit
    type(oscpar_t), allocatable, intent(out) :: op(:)

    integer                              :: n_op
    integer                              :: i
    character(LEN(op%rotation_method))   :: rotation_method
    character(LEN(op%variables_type))    :: variables_type
    character(LEN(op%inner_bound_type))  :: inner_bound_type
    character(LEN(op%outer_bound_type))  :: outer_bound_type
    character(LEN(op%inertia_norm_type)) :: inertia_norm_type
    character(LEN(op%tag_list))          :: tag_list
    logical                              :: reduce_order
    logical                              :: nonadiabatic
    real(WP)                             :: x_ref

    namelist /osc/ x_ref, rotation_method, inner_bound_type, outer_bound_type, variables_type, &
         inertia_norm_type, tag_list, reduce_order, nonadiabatic

    ! Count the number of osc namelists

    rewind(unit)

    n_op = 0

    count_loop : do
       read(unit, NML=osc, END=100)
       n_op = n_op + 1
    end do count_loop

100 continue

    ! Read oscillation parameters

    rewind(unit)

    allocate(op(n_op))

    read_loop : do i = 1,n_op

       x_ref = HUGE(0._WP)

       rotation_method = 'NULL'
       variables_type = 'DZIEM'
       inner_bound_type = 'REGULAR'
       outer_bound_type = 'ZERO'
       inertia_norm_type = 'BOTH'
       tag_list = ''

       nonadiabatic = .FALSE.
       reduce_order = .TRUE.

       read(unit, NML=osc)

       ! Initialize the oscpar

       op(i) = oscpar_t(x_ref=x_ref, &
                        rotation_method=rotation_method, &
                        variables_type=variables_type, &
                        inner_bound_type=inner_bound_type, &
                        outer_bound_type=outer_bound_type, &
                        inertia_norm_type=inertia_norm_type, &
                        tag_list=tag_list, &
                        nonadiabatic=nonadiabatic, &
                        reduce_order=reduce_order)

    end do read_loop

    ! Finish

    return

  end subroutine read_oscpar

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (op, root_rank)

    type(oscpar_t), intent(inout) :: op$ARRAY_SPEC($RANK)
    integer, intent(in)           :: root_rank

    ! Broadcast the oscpar_t

    call bcast(op%x_ref, root_rank)

    call bcast(op%rotation_method, root_rank)
    call bcast(op%variables, root_rank)
    call bcast(op%inner_bound, root_rank)
    call bcast(op%outer_bound, root_rank)
    call bcast(op%inertia_norm, root_rank)
    call bcast(op%tag_list, root_rank)

    call bcast(np%nonadiabatic, root_rank)
    call bcast(np%reduce_order, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(oscpar_t),0)
  $BCAST_ALLOC(type(oscpar_t),1)

  $endif

end module gyre_oscpar
