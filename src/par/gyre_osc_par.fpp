! Module   : gyre_osc_par
! Purpose  : oscillation parameters
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
$include 'core_parallel.inc'

module gyre_osc_par

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: osc_par_t
     real(WP)                :: x_ref = 1._WP
     real(WP)                :: x_atm = -1._WP
     real(WP)                :: alpha_grv = 1._WP
     real(WP)                :: alpha_thm = 1._WP
     real(WP)                :: alpha_hfl = 1._WP
     real(WP)                :: alpha_gam = 1._WP
     real(WP)                :: alpha_pi = 1._WP
     real(WP)                :: alpha_kap = 1._WP
     character(64)           :: variables_set = 'GYRE'
     character(64)           :: inner_bound = 'REGULAR'
     character(64)           :: outer_bound = 'VACUUM'
     character(64)           :: outer_bound_cutoff = ''
     character(64)           :: outer_bound_branch = 'E_NEG'
     character(64)           :: inertia_norm = 'BOTH'
     character(64)           :: time_factor = 'OSC'
     character(64)           :: conv_scheme = 'FROZEN_PESNELL_1'
     character(64)           :: zeta_scheme = 'PESNELL'
     character(64)           :: deps_source = 'MODEL'
     character(FILENAME_LEN) :: deps_file = ''
     character(256)          :: deps_file_format = ''
     character(2048)         :: tag_list = ''
     logical                 :: adiabatic = .TRUE.
     logical                 :: nonadiabatic = .FALSE.
     logical                 :: quasiad_eigfuncs = .FALSE.
     logical                 :: eddington_approx = .FALSE.
     logical                 :: reduce_order = .TRUE.
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

    integer                               :: n_os_p
    integer                               :: i
    real(WP)                              :: x_ref
    real(WP)                              :: x_atm
    real(WP)                              :: alpha_grv
    real(WP)                              :: alpha_thm
    real(WP)                              :: alpha_hfl
    real(WP)                              :: alpha_gam
    real(WP)                              :: alpha_pi
    real(WP)                              :: alpha_kap
    character(LEN(os_p%variables_set))    :: variables_set
    character(LEN(os_p%inner_bound))      :: inner_bound
    character(LEN(os_p%outer_bound))      :: outer_bound
    character(LEN(os_p%outer_bound))      :: outer_bound_cutoff 
    character(LEN(os_p%outer_bound))      :: outer_bound_branch
    character(LEN(os_p%inertia_norm))     :: inertia_norm
    character(LEN(os_p%time_factor))      :: time_factor
    character(LEN(os_p%conv_scheme))      :: conv_scheme
    character(LEN(os_p%zeta_scheme))      :: zeta_scheme
    character(LEN(os_p%deps_source))      :: deps_source
    character(LEN(os_p%deps_file))        :: deps_file
    character(LEN(os_p%deps_file_format)) :: deps_file_format
    character(LEN(os_p%tag_list))         :: tag_list
    logical                               :: adiabatic
    logical                               :: nonadiabatic
    logical                               :: quasiad_eigfuncs
    logical                               :: eddington_approx
    logical                               :: reduce_order

    namelist /osc/ x_ref, x_atm, alpha_grv, alpha_thm, alpha_hfl, &
         alpha_gam, alpha_pi, alpha_kap, &
         inner_bound, outer_bound, outer_bound_cutoff, outer_bound_branch, &
         variables_set, inertia_norm, time_factor, &
         conv_scheme, zeta_scheme, deps_source, deps_file, deps_file_format, &
         tag_list, adiabatic, nonadiabatic, quasiad_eigfuncs, &
         eddington_approx, reduce_order

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

       ! Set default values

       os_p(i) = osc_par_t()

       x_ref = os_p(i)%x_ref
       x_atm = os_p(i)%x_atm
       alpha_grv = os_p(i)%alpha_grv
       alpha_thm = os_p(i)%alpha_thm
       alpha_hfl = os_p(i)%alpha_hfl
       alpha_gam = os_p(i)%alpha_gam
       alpha_pi = os_p(i)%alpha_pi
       alpha_kap = os_p(i)%alpha_kap
       variables_set = os_p(i)%variables_set
       inner_bound = os_p(i)%inner_bound
       outer_bound = os_p(i)%outer_bound
       outer_bound_cutoff = os_p(i)%outer_bound_cutoff
       outer_bound_branch = os_p(i)%outer_bound_branch
       inertia_norm = os_p(i)%inertia_norm
       time_factor = os_p(i)%time_factor
       conv_scheme = os_p(i)%conv_scheme
       zeta_scheme = os_p(i)%zeta_scheme
       deps_source = os_p(i)%deps_source
       deps_file = os_p(i)%deps_file
       deps_file_format = os_p(i)%deps_file_format
       tag_list = os_p(i)%tag_list
       adiabatic = os_p(i)%adiabatic
       nonadiabatic = os_p(i)%nonadiabatic
       quasiad_eigfuncs = os_p(i)%quasiad_eigfuncs
       eddington_approx = os_p(i)%eddington_approx
       reduce_order = os_p(i)%reduce_order

       ! Read the namelist

       read(unit, NML=osc)

       ! Store read values

       os_p(i)%x_ref = x_ref
       os_p(i)%x_atm = x_atm
       os_p(i)%alpha_grv = alpha_grv
       os_p(i)%alpha_thm = alpha_thm
       os_p(i)%alpha_hfl = alpha_hfl
       os_p(i)%alpha_gam = alpha_gam
       os_p(i)%alpha_pi = alpha_pi
       os_p(i)%alpha_kap = alpha_kap
       os_p(i)%variables_set = variables_set
       os_p(i)%inner_bound = inner_bound
       os_p(i)%outer_bound = outer_bound
       os_p(i)%outer_bound_cutoff = outer_bound_cutoff
       os_p(i)%outer_bound_branch = outer_bound_branch
       os_p(i)%inertia_norm = inertia_norm
       os_p(i)%time_factor = time_factor
       os_p(i)%conv_scheme = conv_scheme
       os_p(i)%zeta_scheme = zeta_scheme
       os_p(i)%deps_source = deps_source
       os_p(i)%deps_file = deps_file
       os_p(i)%deps_file_format = deps_file_format
       os_p(i)%tag_list = tag_list
       os_p(i)%adiabatic = adiabatic
       os_p(i)%nonadiabatic = nonadiabatic
       os_p(i)%quasiad_eigfuncs = quasiad_eigfuncs
       os_p(i)%eddington_approx = eddington_approx
       os_p(i)%reduce_order = reduce_order

    end do read_loop

    ! Finish

    return

  end subroutine read_osc_par

  !****

  $if ($MPI)

  subroutine bcast_0_ (os_p, root_rank)

    type(osc_par_t), intent(inout) :: os_p
    integer, intent(in)            :: root_rank

    ! Broadcast the osc_par_t

    call bcast(os_p%x_ref, root_rank)
    call bcast(os_p%x_atm, root_rank)
    call bcast(os_p%alpha_grv, root_rank)
    call bcast(os_p%alpha_thm, root_rank)
    call bcast(os_p%alpha_hfl, root_rank)
    call bcast(os_p%alpha_gam, root_rank)
    call bcast(os_p%alpha_pi, root_rank)
    call bcast(os_p%alpha_kap, root_rank)

    call bcast(os_p%variables_set, root_rank)
    call bcast(os_p%inner_bound, root_rank)
    call bcast(os_p%outer_bound, root_rank)
    call bcast(os_p%outer_bound_cutoff, root_rank)
    call bcast(os_p%outer_bound_branch, root_rank)
    call bcast(os_p%inertia_norm, root_rank)
    call bcast(os_p%time_factor, root_rank)
    call bcast(os_p%conv_scheme, root_rank)
    call bcast(os_p%deps_source, root_rank)
    call bcast(os_p%deps_file, root_rank)
    call bcast(os_p%deps_file_format, root_rank)
    call bcast(os_p%tag_list, root_rank)

    call bcast(os_p%adiabatic, root_rank)
    call bcast(os_p%nonadiabatic, root_rank)
    call bcast(os_p%quasiad_eigfuncs, root_rank)
    call bcast(os_p%eddington_approx, root_rank)
    call bcast(os_p%reduce_order, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(osc_par_t),1)

  $BCAST_ALLOC(type(osc_par_t),0)
  $BCAST_ALLOC(type(osc_par_t),1)

  $endif

end module gyre_osc_par
