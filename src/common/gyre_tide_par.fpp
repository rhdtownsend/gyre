! Module   : gyre_tide_par
! Purpose  : tidal parameters
!
! Copyright 2018 Rich Townsend
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

module gyre_tide_par

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: tide_par_t
     real(WP) :: R_a
     real(WP) :: q
     real(WP) :: e
     real(WP) :: omega_static
     integer  :: l_max
     integer  :: k_max
  end type tide_par_t

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

  public :: tide_par_t
  public :: read_tide_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  subroutine read_tide_par (unit, td_p)

    integer, intent(in)                        :: unit
    type(tide_par_t), allocatable, intent(out) :: td_p(:)

    integer  :: n_td_p
    integer  :: i
    real(WP) :: R_a
    real(WP) :: q
    real(WP) :: e
    real(WP) :: omega_static
    integer  :: l_max
    integer  :: k_max

    namelist /tide/ R_a, q, e, omega_static, l_max, k_max

    ! Count the number of tide namelists

    rewind(unit)

    n_td_p = 0

    count_loop : do
       read(unit, NML=tide, END=100)
       n_td_p = n_td_p + 1
    end do count_loop

100 continue

    ! Read tide parameters

    rewind(unit)

    allocate(td_p(n_td_p))

    read_loop : do i = 1,n_td_p

       R_a = 0.2_WP
       q = 1._WP
       e = 0.5_WP
       omega_static = 0._WP

       l_max = 4
       k_max = 20

       read(unit, NML=tide)

       ! Initialize the tide_par

       td_p(i) = tide_par_t(R_a=R_a, &
                            q=q, &
                            e=e, &
                            omega_static=omega_static, &
                            l_max=l_max, &
                            k_max=k_max)

    end do read_loop

    ! Finish

    return

  end subroutine read_tide_par

  !****

  $if ($MPI)

  subroutine bcast_0_ (td_p, root_rank)

    type(tide_par_t), intent(inout) :: td_p
    integer, intent(in)             :: root_rank

    ! Broadcast the tide_par_t

    call bcast(td_p%R_a, root_rank)
    call bcast(td_p%q, root_rank)
    call bcast(td_p%e, root_rank)
    call bcast(td_p%Omega_static, root_rank)

    call bcast(td_p%l_max, root_rank)
    call bcast(td_p%k_max, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(tide_par_t),1)

  $BCAST_ALLOC(type(tide_par_t),0)
  $BCAST_ALLOC(type(tide_par_t),1)

  $endif

end module gyre_tide_par
