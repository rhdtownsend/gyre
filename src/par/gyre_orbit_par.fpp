! Module   : gyre_orbit_par
! Purpose  : orbit parameters
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

module gyre_orbit_par

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: orbit_par_t
     real(WP)        :: Omega_orb = 1._WP
     real(WP)        :: q = 1._WP
     real(WP)        :: e = 0.5_WP
     real(WP)        :: t_0 = 0._WP
     real(WP)        :: sync_fraction = 0._WP
     logical         :: sync_rot = .FALSE.
     character(64)   :: Omega_orb_units = 'NONE'
     character(2048) :: tag_list = ''
  end type orbit_par_t

 ! Access specifiers

  private

  public :: orbit_par_t
  public :: read_orbit_par

  ! Procedures

contains

  subroutine read_orbit_par (unit, or_p)

    integer, intent(in)                         :: unit
    type(orbit_par_t), allocatable, intent(out) :: or_p(:)

    integer                              :: n_or_p
    integer                              :: i
    real(WP)                             :: Omega_orb
    real(WP)                             :: q
    real(WP)                             :: e
    real(WP)                             :: t_0
    real(WP)                             :: sync_fraction
    logical                              :: sync_rot
    character(LEN(or_p%Omega_orb_units)) :: Omega_orb_units
    character(LEN(or_p%tag_list))        :: tag_list

    namelist /orbit/ Omega_orb, q, e, t_0, sync_fraction, sync_rot, &
         Omega_orb_units, tag_list

    ! Count the number of orbit namelists

    rewind(unit)

    n_or_p = 0

    count_loop : do
       read(unit, NML=orbit, END=100)
       n_or_p = n_or_p + 1
    end do count_loop

100 continue

    ! Read orbit parameters

    rewind(unit)

    allocate(or_p(n_or_p))

    read_loop : do i = 1,n_or_p

       ! Set default values

       or_p(i) = orbit_par_t()

       Omega_orb = or_p(i)%Omega_orb
       q = or_p(i)%q
       e = or_p(i)%e
       t_0 = or_p(i)%t_0
       sync_fraction = or_p(i)%sync_fraction

       sync_rot = or_p(i)%sync_rot

       Omega_orb_units = or_p(i)%Omega_orb_units
       tag_list = or_p(i)%tag_list

       ! Read the namelist

       read(unit, NML=orbit)

       ! Store read values

       or_p(i)%Omega_orb = Omega_orb
       or_p(i)%q = q
       or_p(i)%e = e
       or_p(i)%t_0 = t_0
       or_p(i)%sync_fraction = sync_fraction

       or_p(i)%sync_rot = sync_rot

       or_p(i)%Omega_orb_units = Omega_orb_units
       or_p(i)%tag_list = tag_list

    end do read_loop

    ! Finish

    return

  end subroutine read_orbit_par

end module gyre_orbit_par
