! Program  : gyre_orbit
! Purpose  : secular orbital evolution code
!
! Copyright 2018-2019 Rich Townsend & The GYRE Team
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

program gyre_orbit

  ! Uses

  use core_kinds, only : WP
  use core_hgroup
  use core_parallel
  use core_system

  use gyre_context
  use gyre_constants
  use gyre_evol_model
  use gyre_func
  use gyre_grid_par
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_scan
  use gyre_scan_par
  use gyre_search
  use gyre_tide
  use gyre_tide_par
  use gyre_tide_util
  use gyre_util
  use gyre_version
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable      :: filename
  integer                        :: unit
  type(model_par_t)              :: ml_p
  type(osc_par_t), allocatable   :: os_p(:)
  type(num_par_t), allocatable   :: nm_p(:)
  type(grid_par_t), allocatable  :: gr_p(:)
  type(scan_par_t), allocatable  :: sc_p(:)
  type(tide_par_t), allocatable  :: td_p(:)
  class(model_t), pointer        :: ml => null()
  type(mode_par_t)               :: md_p
  type(context_t)                :: cx
  real(WP), allocatable          :: Omega_orb(:)
  integer                        :: n_Omega_orb
  real(WP), allocatable          :: a_dot(:)
  real(WP), allocatable          :: e_dot(:)
  real(WP), allocatable          :: o_dot(:)
  real(WP), allocatable          :: J_dot(:)
  integer                        :: i
  type(hgroup_t)                 :: hg

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_orbit <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_orbit ['//VERSION//']', '-')
100  format(A)

     if (check_log_level('DEBUG')) then
        write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
        write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110     format(A,1X,A)
     endif

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

     write(OUTPUT_UNIT, *)

  endif

  ! Read the namelist file

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_scan_par(unit, sc_p)
  call read_tide_par(unit, td_p)

  close(unit)

  $ASSERT(SIZE(os_p) == 1,Must be exactly one osc parameter)
  $ASSERT(SIZE(nm_p) == 1,Must be exactly one num parameter)
  $ASSERT(SIZE(gr_p) == 1,Must be exactly one grid parameter)
  $ASSERT(SIZE(sc_p) >= 1,Must be at least one scan parameter)
  $ASSERT(SIZE(td_p) == 1,Must be exactly one tide parameter)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '-')
  endif

  ml => model_t(ml_p)

  ! Create the reference mode_par_t for setting up contexts, frequency
  ! searches

  md_p = mode_par_t(l=td_p(1)%l_ref, m=td_p(1)%m_ref)

  ! Set up the context

  cx = context_t(ml, gr_p(1), md_p, os_p(1))

  ! Set up the frequency array

  call build_scan(cx, md_p, os_p(1), sc_p, Omega_orb)

  ! Allocate secular rate-of-change arrays

  n_Omega_orb = SIZE(Omega_orb)

  allocate(a_dot(n_Omega_orb))
  allocate(e_dot(n_Omega_orb))
  allocate(o_dot(n_Omega_orb))
  allocate(J_dot(n_Omega_orb))

  ! Loop over orbital frequencies

  Omega_orb_loop : do i = 1, n_Omega_orb

     ! Initialize the secular arrays

     a_dot(i) = 0._WP
     e_dot(i) = 0._WP
     o_dot(i) = 0._WP
     J_dot(i) = 0._WP

     ! Add in contributions from each tidal component

     call eval_tide(ml, process_wave, Omega_orb(i), os_p(1), nm_p(1), gr_p(1), td_p(1))

  end do Omega_orb_loop

  ! Write out results

  hg = hgroup_t('orbit.h5', CREATE_FILE)

  select type (ml)
  class is (evol_model_t)
     call write_attr(hg, 'M_star', ml%M_star)
     call write_attr(hg, 'R_star', ml%R_star)
     call write_attr(hg, 'L_star', ml%L_star)
  end select

  call write_attr(hg, 'q', td_p(1)%q)
  call write_attr(hg, 'e', td_p(1)%e)
  call write_attr(hg, 'l_max', td_p(1)%l_max)
  call write_attr(hg, 'k_max', td_p(1)%k_max)

  call write_dset(hg, 'Omega_orb', Omega_orb)

  call write_dset(hg, 'a_dot', a_dot)
  call write_dset(hg, 'e_dot', e_dot)
  call write_dset(hg, 'o_dot', o_dot)
  call write_dset(hg, 'J_dot', J_dot)

  call hg%final()

  ! Clean up

  deallocate(ml)

  ! Finish

  call final_parallel()

contains

  subroutine process_wave (wv, k)

    type(wave_t), intent(in) :: wv
    integer, intent(in)      :: k

    real(WP)    :: e
    real(WP)    :: q
    real(WP)    :: R_a
    real(WP)    :: eps_tide
    integer     :: l
    integer     :: m
    real(WP)    :: c
    real(WP)    :: kappa
    complex(WP) :: F
    real(WP)    :: gamma

    ! Set up orbital parameters

    e = td_p(1)%e
    q = td_p(1)%q

    R_a = (Omega_orb(i)**2/(1._WP + td_p(1)%q))**(1._WP/3._WP)
    
    eps_tide = (R_a)**3*q

    ! Evaluate the summation coefficient

    l = wv%l
    m = wv%m

    if (k == 0) then

       if (m == 0) then
          kappa = 0.5_WP
       elseif (m >= 1) then
          kappa = 1._WP
       else
          kappa = 0._WP
       endif

    else

       kappa = 1._WP

    endif

    ! Evaluate the tidal potential coefficient

    c = tidal_c(R_a, e, l, m, k)

    print *,'R_a, l, m, k:',R_a, e, l, m, k

    if (kappa /= 0._WP .AND. c /= 0._WP) then

       ! Evaluate the response function
 
       F = -0.5_WP*(SQRT(4._WP*PI)*wv%eul_phi(wv%n_k)/(eps_tide*c) + 1._WP)

       ! Evaluate the summation weight

       ! Accumulate the (dimensionless) secular rates-of-change
       ! elements

       gamma = ATAN2(AIMAG(F), REAL(F))

       a_dot(i) = a_dot(i) + 4._WP*Omega_orb(i)*(q/R_a)*(R_a)**(l+3)* &
            kappa*ABS(F)*SIN(gamma)*secular_G_2(R_a, e, l, m, k)

       e_dot(i) = e_dot(i) + 4._WP*Omega_orb(i)*q*(R_a)**(l+3)* &
            kappa*ABS(F)*SIN(gamma)*secular_G_3(R_a, e, l, m, k)

       o_dot(i) = o_dot(i) + 4._WP*Omega_orb(i)*q*(R_a)**(l+3)* &
            kappa*ABS(F)*COS(gamma)*secular_G_1(R_a, e, l, m, k)

       ! Accumulate the torque

       J_dot(i) = J_dot(i) + wv%tau_ss()

    endif

    ! Finish

    return

  end subroutine process_wave

end program gyre_orbit
