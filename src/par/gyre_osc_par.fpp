! Module   : gyre_osc_par
! Purpose  : oscillation parameters
!
! Copyright 2013-2018 Rich Townsend
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
  use core_constants, only : FILENAME_LEN
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: osc_par_t
     real(WP)                :: x_ref = 1._WP
     real(WP)                :: alpha_th = 1._WP
     character(64)           :: rotation_method = 'DOPPLER'
     character(64)           :: variables_set = 'GYRE'
     character(64)           :: inner_bound = 'REGULAR'
     character(64)           :: outer_bound = 'VACUUM'
     character(64)           :: outer_bound_for_cutoff = 'UNNO'
     character(64)           :: inertia_norm = 'BOTH'
     character(64)           :: time_factor = 'OSC'
     character(64)           :: conv_scheme = 'FROZEN_PESNELL_1'
     character(64)           :: int_scheme = 'PESNELL'
     character(64)           :: deps_scheme = 'MODEL'
     character(FILENAME_LEN) :: deps_file = ''
     character(256)          :: deps_file_format = ''
     character(2048)         :: tag_list = ''
     logical                 :: adiabatic = .TRUE.
     logical                 :: nonadiabatic = .FALSE.
     logical                 :: quasiad_eigfuncs = .FALSE.
     logical                 :: cowling_approx = .FALSE.
     logical                 :: nar_approx = .FALSE.
     logical                 :: narf_approx = .FALSE.
     logical                 :: eddington_approx = .FALSE.
     logical                 :: complex_lambda = .FALSE.
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
    real(WP)                              :: alpha_th
    character(LEN(os_p%rotation_method))  :: rotation_method
    character(LEN(os_p%variables_set))    :: variables_set
    character(LEN(os_p%inner_bound))      :: inner_bound
    character(LEN(os_p%outer_bound))      :: outer_bound
    character(LEN(os_p%outer_bound))      :: outer_bound_for_cutoff
    character(LEN(os_p%inertia_norm))     :: inertia_norm
    character(LEN(os_p%time_factor))      :: time_factor
    character(LEN(os_p%conv_scheme))      :: conv_scheme
    character(LEN(os_p%int_scheme))       :: int_scheme
    character(LEN(os_p%deps_scheme))      :: deps_scheme
    character(LEN(os_p%deps_file))        :: deps_file
    character(LEN(os_p%deps_file_format)) :: deps_file_format
    character(LEN(os_p%tag_list))         :: tag_list
    logical                               :: adiabatic
    logical                               :: nonadiabatic
    logical                               :: quasiad_eigfuncs
    logical                               :: cowling_approx
    logical                               :: nar_approx
    logical                               :: narf_approx
    logical                               :: eddington_approx
    logical                               :: complex_lambda
    logical                               :: reduce_order

    namelist /osc/ x_ref, alpha_th, rotation_method, inner_bound, outer_bound, &
         outer_bound_for_cutoff, variables_set, inertia_norm, time_factor, &
         conv_scheme, int_scheme, deps_scheme, deps_file, deps_file_format, &
         tag_list, adiabatic, nonadiabatic, quasiad_eigfuncs, &
         cowling_approx, nar_approx, narf_approx, eddington_approx, &
         complex_lambda, reduce_order

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
       alpha_th = os_p(i)%alpha_th
       rotation_method = os_p(i)%rotation_method
       variables_set = os_p(i)%variables_set
       inner_bound = os_p(i)%inner_bound
       outer_bound = os_p(i)%outer_bound
       outer_bound_for_cutoff = os_p(i)%outer_bound_for_cutoff
       inertia_norm = os_p(i)%inertia_norm
       time_factor = os_p(i)%time_factor
       conv_scheme = os_p(i)%conv_scheme
       int_scheme = os_p(i)%int_scheme
       deps_scheme = os_p(i)%deps_scheme
       deps_file = os_p(i)%deps_file
       deps_file_format = os_p(i)%deps_file_format
       tag_list = os_p(i)%tag_list
       adiabatic = os_p(i)%adiabatic
       nonadiabatic = os_p(i)%nonadiabatic
       quasiad_eigfuncs = os_p(i)%quasiad_eigfuncs
       cowling_approx = os_p(i)%cowling_approx
       nar_approx = os_p(i)%nar_approx
       narf_approx = os_p(i)%narf_approx
       eddington_approx = os_p(i)%eddington_approx
       complex_lambda = os_p(i)%complex_lambda
       reduce_order = os_p(i)%reduce_order

       ! Read the namelist

       read(unit, NML=osc)

       ! Store read values

       os_p(i)%x_ref = x_ref
       os_p(i)%alpha_th = alpha_th
       os_p(i)%rotation_method = rotation_method
       os_p(i)%variables_set = variables_set
       os_p(i)%inner_bound = inner_bound
       os_p(i)%outer_bound = outer_bound
       os_p(i)%outer_bound_for_cutoff = outer_bound_for_cutoff
       os_p(i)%inertia_norm = inertia_norm
       os_p(i)%time_factor = time_factor
       os_p(i)%conv_scheme = conv_scheme
       os_p(i)%int_scheme = int_scheme
       os_p(i)%deps_scheme = deps_scheme
       os_p(i)%deps_file = deps_file
       os_p(i)%deps_file_format = deps_file_format
       os_p(i)%tag_list = tag_list
       os_p(i)%adiabatic = adiabatic
       os_p(i)%nonadiabatic = nonadiabatic
       os_p(i)%quasiad_eigfuncs = quasiad_eigfuncs
       os_p(i)%cowling_approx = cowling_approx
       os_p(i)%nar_approx = nar_approx
       os_p(i)%narf_approx = narf_approx
       os_p(i)%eddington_approx = eddington_approx
       os_p(i)%complex_lambda = complex_lambda
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
    call bcast(os_p%alpha_th, root_rank)

    call bcast(os_p%rotation_method, root_rank)
    call bcast(os_p%variables_set, root_rank)
    call bcast(os_p%inner_bound, root_rank)
    call bcast(os_p%outer_bound, root_rank)
    call bcast(os_p%outer_bound_for_cutoff, root_rank)
    call bcast(os_p%inertia_norm, root_rank)
    call bcast(os_p%time_factor, root_rank)
    call bcast(os_p%conv_scheme, root_rank)
    call bcast(os_p%deps_scheme, root_rank)
    call bcast(os_p%deps_file, root_rank)
    call bcast(os_p%deps_file_format, root_rank)
    call bcast(os_p%tag_list, root_rank)

    call bcast(os_p%adiabatic, root_rank)
    call bcast(os_p%nonadiabatic, root_rank)
    call bcast(os_p%quasiad_eigfuncs, root_rank)
    call bcast(os_p%cowling_approx, root_rank)
    call bcast(os_p%nar_approx, root_rank)
    call bcast(os_p%narf_approx, root_rank)
    call bcast(os_p%eddington_approx, root_rank)
    call bcast(os_p%complex_lambda, root_rank)
    call bcast(os_p%reduce_order, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(osc_par_t),1)

  $BCAST_ALLOC(type(osc_par_t),0)
  $BCAST_ALLOC(type(osc_par_t),1)

  $endif

end module gyre_osc_par
