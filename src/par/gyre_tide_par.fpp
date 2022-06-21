! Module   : gyre_tide_par
! Purpose  : tidal parameters
!
! Copyright 2018-2022 Rich Townsend & The GYRE Team
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

module gyre_tide_par

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: tide_par_t
     real(WP)      :: Psi_o_thresh = 0._WP
     real(WP)      :: omega_c_thresh = 0._WP
     integer       :: l_min = 2
     integer       :: l_max = 4
     integer       :: m_min = -HUGE(0)
     integer       :: m_max = HUGE(0)
     integer       :: k_min = -10
     integer       :: k_max = 10
     character(64) :: tag = ''
  end type tide_par_t

 ! Access specifiers

  private

  public :: tide_par_t
  public :: read_tide_par

  ! Procedures

contains

  subroutine read_tide_par (unit, td_p)

    integer, intent(in)                        :: unit
    type(tide_par_t), allocatable, intent(out) :: td_p(:)

    integer  :: n_td_p
    integer  :: i
    real(WP) :: Psi_o_thresh
    real(WP) :: omega_c_thresh
    integer  :: l_min
    integer  :: l_max
    integer  :: m_min
    integer  :: m_max
    integer  :: k_min
    integer  :: k_max

    namelist /tide/ Psi_o_thresh, omega_c_thresh, &
         l_min, l_max, m_min, m_max, k_min, k_max

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

       Psi_o_thresh = td_p(i)%Psi_o_thresh
       omega_c_thresh = td_p(i)%omega_c_thresh

       l_min = td_p(i)%l_min
       l_max = td_p(i)%l_max
       m_min = td_p(i)%m_min
       m_max = td_p(i)%m_max
       k_min = td_p(i)%k_min
       k_max = td_p(i)%k_max

       ! Read the namelist

       read(unit, NML=tide)

       ! Store read values

       td_p(i)%Psi_o_thresh = Psi_o_thresh
       td_p(i)%omega_c_thresh = omega_c_thresh

       td_p(i)%l_min = l_min
       td_p(i)%l_max = l_max
       td_p(i)%m_min = m_min
       td_p(i)%m_max = m_max
       td_p(i)%k_min = k_min
       td_p(i)%k_max = k_max

    end do read_loop

    ! Finish

    return

  end subroutine read_tide_par

end module gyre_tide_par
