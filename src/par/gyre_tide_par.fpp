! Module   : gyre_tide_par
! Purpose  : tidal parameters
!
! Copyright 2018-2020 Rich Townsend & The GYRE Team
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
     real(WP) :: omega_static = 0._WP
     integer  :: l_ref = 0
     integer  :: m_ref = 0
     integer  :: l_max = 4
     integer  :: m_min = -HUGE(0)
     integer  :: m_max = HUGE(0)
     integer  :: k_min = -20
     integer  :: k_max = 20
     logical  :: combine_k = .TRUE.
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
    real(WP) :: omega_static
    integer  :: l_ref
    integer  :: m_ref
    integer  :: l_max
    integer  :: m_min
    integer  :: m_max
    integer  :: k_min
    integer  :: k_max
    logical  :: combine_k

    namelist /tide/ omega_static, l_ref, m_ref, l_max, m_min, m_max, &
         k_min, k_max, combine_k

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

       ! Set default values

       td_p(i) = tide_par_t()

       omega_static = td_p(i)%omega_static

       l_ref = td_p(i)%l_ref
       m_ref = td_p(i)%m_ref
       l_max = td_p(i)%l_max
       m_min = td_p(i)%m_min
       m_max = td_p(i)%m_max
       k_min = td_p(i)%k_min
       k_max = td_p(i)%k_max

       combine_k = td_p(i)%combine_k

       ! Read the namelist

       read(unit, NML=tide)

       ! Store read values

       td_p(i)%omega_static = omega_static

       td_p(i)%l_ref = l_ref
       td_p(i)%m_ref = m_ref
       td_p(i)%l_max = l_max
       td_p(i)%m_min = m_min
       td_p(i)%m_max = m_max
       td_p(i)%k_min = k_min
       td_p(i)%k_max = k_max

       td_p(i)%combine_k = combine_k

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

    call bcast(td_p%omega_static, root_rank)

    call bcast(td_p%l_ref, root_rank)
    call bcast(td_p%m_ref, root_rank)
    call bcast(td_p%l_max, root_rank)
    call bcast(td_p%m_min, root_rank)
    call bcast(td_p%m_max, root_rank)
    call bcast(td_p%k_min, root_rank)
    call bcast(td_p%k_max, root_rank)

    call bcast(td_p%combine_k, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(tide_par_t),1)

  $BCAST_ALLOC(type(tide_par_t),0)
  $BCAST_ALLOC(type(tide_par_t),1)

  $endif

end module gyre_tide_par
