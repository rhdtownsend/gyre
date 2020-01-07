! Module   : gyre_rot_par
! Purpose  : rotation parameters
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
$include 'core_parallel.inc'

module gyre_rot_par

  ! Uses

  use core_kinds
  use core_constants, only : FILENAME_LEN
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: rot_par_t
     real(WP)        :: Omega_rot = 0._WP
     character(64)   :: coriolis_method = 'NULL'
     character(64)   :: Omega_rot_source = 'MODEL'
     character(64)   :: Omega_rot_units = 'NONE'
     character(2048) :: tag_list = ''
     logical         :: complex_lambda = .FALSE.
  end type rot_par_t

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

  public :: rot_par_t
  public :: read_rot_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  subroutine read_rot_par (unit, rt_p)

    integer, intent(in)                       :: unit
    type(rot_par_t), allocatable, intent(out) :: rt_p(:)

    integer                               :: n_rt_p
    integer                               :: i
    real(WP)                              :: Omega_rot
    character(LEN(rt_p%coriolis_method))  :: coriolis_method
    character(LEN(rt_p%Omega_rot_source)) :: Omega_rot_source
    character(LEN(rt_p%Omega_rot_units))  :: Omega_rot_units
    character(LEN(rt_p%tag_list))         :: tag_list
    logical                               :: complex_lambda

    namelist /rot/ Omega_rot, coriolis_method, Omega_rot_source, Omega_rot_units, &
         tag_list, complex_lambda

    ! Count the number of rot namelists

    rewind(unit)

    n_rt_p = 0

    count_loop : do
       read(unit, NML=rot, END=100)
       n_rt_p = n_rt_p + 1
    end do count_loop

100 continue

    ! Read rotation parameters

    rewind(unit)

    allocate(rt_p(n_rt_p))

    read_loop : do i = 1,n_rt_p

       ! Set default values

       rt_p(i) = rot_par_t()

       Omega_rot = rt_p(i)%Omega_rot
       coriolis_method = rt_p(i)%coriolis_method
       Omega_rot_source = rt_p(i)%Omega_rot_source
       Omega_rot_units = rt_p(i)%Omega_rot_units
       tag_list = rt_p(i)%tag_list
       complex_lambda = rt_p(i)%complex_lambda

       ! Read the namelist

       read(unit, NML=rot)

       ! Store read values

       rt_p(i)%Omega_rot = Omega_rot
       rt_p(i)%coriolis_method = coriolis_method
       rt_p(i)%Omega_rot_source = Omega_rot_source
       rt_p(i)%Omega_rot_units = Omega_rot_units
       rt_p(i)%tag_list = tag_list
       rt_p(i)%complex_lambda = complex_lambda

    end do read_loop

    ! Finish

    return

  end subroutine read_rot_par

  !****

  $if ($MPI)

  subroutine bcast_0_ (rt_p, root_rank)

    type(rot_par_t), intent(inout) :: rt_p
    integer, intent(in)            :: root_rank

    ! Broadcast the rot_par_t

    call bcast(rt_p%Omega_rot, root_rank)

    call bcast(rt_p%coriolis_method, root_rank)
    call bcast(rt_p%Omega_rot_source, root_rank)
    call bcast(rt_p%Omega_rot_units, root_rank)
    call bcast(rt_p%tag_list, root_rank)

    call bcast(rt_p%complex_lambda, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(rot_par_t),1)

  $BCAST_ALLOC(type(rot_par_t),0)
  $BCAST_ALLOC(type(rot_par_t),1)

  $endif

end module gyre_rot_par
