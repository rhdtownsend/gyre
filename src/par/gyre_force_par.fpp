! Module   : gyre_force_par
! Purpose  : forcing parameters
!
! Copyright 2019 Rich Townsend
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

module gyre_force_par

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: force_par_t
     real(WP)        :: Phi = 0._WP
     real(WP)        :: q = 1._WP
     real(WP)        :: e = 0._WP
     integer         :: k = 1
     character(64)   :: force_type = 'FIXED'
     character(2048) :: tag_list = ''
  end type force_par_t

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

  public :: force_par_t
  public :: read_force_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

!****

  subroutine read_force_par (unit, fr_p)

    integer, intent(in)                         :: unit
    type(force_par_t), allocatable, intent(out) :: fr_p(:)

    integer                         :: n_fr_p
    integer                         :: i
    real(WP)                        :: Phi
    real(WP)                        :: q
    real(WP)                        :: e
    integer                         :: k
    character(LEN(fr_p%force_type)) :: force_type
    character(LEN(fr_p%tag_list))   :: tag_list

    namelist /force/ Phi, q, e, k, force_type, tag_list

    ! Count the number of force namelists

    rewind(unit)

    n_fr_p = 0

    count_loop : do
       read(unit, NML=force, END=100)
       n_fr_p = n_fr_p + 1
    end do count_loop

100 continue

    ! Read force parameters

    rewind(unit)

    allocate(fr_p(n_fr_p))

    read_loop : do i = 1,n_fr_p

       ! Set default values

       fr_p(i) = force_par_t()

       Phi = fr_p(i)%Phi
       q = fr_p(i)%q
       e = fr_p(i)%e
       k = fr_p(i)%k
       force_type = fr_p(i)%force_type
       tag_list = fr_p(i)%tag_list

       ! Read the namelist

       read(unit, NML=force)

       ! Store read values

       fr_p(i)%Phi = Phi
       fr_p(i)%q = q
       fr_p(i)%e = e
       fr_p(i)%k = k
       fr_p(i)%force_type = force_type
       fr_p(i)%tag_list = tag_list

    end do read_loop

    ! Finish

    return

  end subroutine read_force_par

  !****

  $if ($MPI)

  subroutine bcast_0_ (fr_p, root_rank)

    type(force_par_t), intent(inout) :: fr_p
    integer, intent(in)              :: root_rank

    ! Broadcast the force_par_t

    call bcast(fr_p%Phi, root_rank)

    call bcast(fr_p%q, root_rank)
    call bcast(fr_p%e, root_rank)
    call bcast(fr_p%k, root_rank)

    call bcast(fr_p%force_type, root_rank)

    call bcast(fr_p%tag_list, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(force_par_t),1)

  $BCAST_ALLOC(type(force_par_t),0)
  $BCAST_ALLOC(type(force_par_t),1)

  $endif

end module gyre_force_par
