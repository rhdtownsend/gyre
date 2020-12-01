! Program  : gyre_orbit
! Purpose  : secular orbital evolution code
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

program gyre_orbit

  ! Uses

  use core_kinds, only : WP
  use core_hgroup
  use core_parallel
  use core_system

  use gyre_constants
  use gyre_evol_model
  use gyre_func
  use gyre_grid_par
  use gyre_math
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_num_par
  use gyre_orbit_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_rot_par
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
  type(rot_par_t), allocatable   :: rt_p(:)
  type(num_par_t), allocatable   :: nm_p(:)
  type(grid_par_t), allocatable  :: gr_p(:)
  type(orbit_par_t), allocatable :: or_p(:)
  type(tide_par_t), allocatable  :: td_p(:)
  type(out_par_t)                :: ot_p
  class(model_t), pointer        :: ml => null()
  integer                        :: n_or_p
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
  call init_math()

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
  call read_rot_par(unit, rt_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_orbit_par(unit, or_p)
  call read_tide_par(unit, td_p)
  call read_out_par(unit, 'tide', ot_p)

  close(unit)

  $ASSERT(SIZE(os_p) == 1,Must be exactly one osc parameter)
  $ASSERT(SIZE(rt_p) == 1,Must be exactly one rot parameter)
  $ASSERT(SIZE(nm_p) == 1,Must be exactly one num parameter)
  $ASSERT(SIZE(gr_p) == 1,Must be exactly one grid parameter)
  $ASSERT(SIZE(or_p) >= 1,Must be at least one orbit parameter)
  $ASSERT(SIZE(td_p) == 1,Must be exactly one tide parameter)

  ! Check that GYRE_DIR is set

  $ASSERT(GYRE_DIR /= '',The GYRE_DIR environment variable is not set)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '-')
  endif

  ml => model_t(ml_p)

  ! Allocate secular rate-of-change arrays

  n_or_p = SIZE(or_p)

  allocate(a_dot(n_or_p))
  allocate(e_dot(n_or_p))
  allocate(o_dot(n_or_p))
  allocate(J_dot(n_or_p))

  ! Loop over orbital parameters

  or_p_loop : do i = 1, n_or_p

     ! Initialize the secular arrays

     a_dot(i) = 0._WP
     e_dot(i) = 0._WP
     o_dot(i) = 0._WP
     J_dot(i) = 0._WP

     ! Add in contributions from each tidal component

     call eval_tide(ml, process_wave, os_p(1), rt_p(1), nm_p(1), gr_p(1), or_p(i), td_p(1))

  end do or_p_loop

  ! Write out results

  hg = hgroup_t(ot_p%summary_file, CREATE_FILE)

  select type (ml)
  class is (evol_model_t)
     call write_attr(hg, 'M_star', ml%M_star)
     call write_attr(hg, 'R_star', ml%R_star)
     call write_attr(hg, 'L_star', ml%L_star)
  end select

  call write_attr(hg, 'l_max', td_p(1)%l_max)
  call write_attr(hg, 'k_max', td_p(1)%k_max)
  call write_attr(hg, 'm_min', td_p(1)%m_min)
  call write_attr(hg, 'm_max', td_p(1)%m_max)

  call write_dset(hg, 'Omega_orb', or_p%Omega_orb)
  call write_dset(hg, 'q', or_p%q)
  call write_dset(hg, 'e', or_p%e)
  call write_dset(hg, 't_0', or_p%t_0)
  call write_dset(hg, 'sync_fraction', or_p%sync_fraction)
  call write_dset(hg, 'Omega_orb_units', or_p%Omega_orb_units)

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

    real(WP)    :: Omega_orb
    real(WP)    :: q
    real(WP)    :: e
    real(WP)    :: R_a
    real(WP)    :: eps_tide
    integer     :: l
    integer     :: m
    real(WP)    :: c
    real(WP)    :: kappa
    complex(WP) :: F
    real(WP)    :: gamma

    ! Set up orbital parameters

    Omega_orb = or_p(i)%Omega_orb
    q = or_p(i)%q
    e = or_p(i)%e

    R_a = (Omega_orb**2/(1._WP + or_p(i)%q))**(1._WP/3._WP)
    
    eps_tide = (R_a)**3*q

    ! Evaluate the summation coefficient

    l = wv%l
    m = wv%m

    if (td_p(1)%combine_k) then

       if (k == 0) then

          if (m == 0) then
             kappa = 0.5_WP
          elseif (m >= 1) then
             kappa = 1._WP
          else
             kappa = 0._WP
          endif
          
       elseif (k > 0) then

          kappa = 1._WP

       else

          kappa = 0._WP

       endif

    else

       kappa = 0.5_WP

    endif

    ! Evaluate the tidal potential coefficient

    c = tidal_c(R_a, e, l, m, k)

    if (kappa /= 0._WP .AND. c /= 0._WP) then

       ! Evaluate the response function
 
       F = -0.5_WP*(sqrt(4._WP*PI)*wv%eul_phi(wv%n_k)/(eps_tide*c) + 1._WP)

       ! Evaluate the summation weight

       ! Accumulate the (dimensionless) secular rates-of-change
       ! elements

       gamma = ATAN2(AIMAG(F), REAL(F))

       a_dot(i) = a_dot(i) + 4._WP*Omega_orb*(q/R_a)*(R_a)**(l+3)* &
            kappa*abs(F)*sin(gamma)*secular_G_2(R_a, e, l, m, k)

       e_dot(i) = e_dot(i) + 4._WP*Omega_orb*q*(R_a)**(l+3)* &
            kappa*abs(F)*sin(gamma)*secular_G_3(R_a, e, l, m, k)

       o_dot(i) = o_dot(i) + 4._WP*Omega_orb*q*(R_a)**(l+3)* &
            kappa*abs(F)*cos(gamma)*secular_G_1(R_a, e, l, m, k)

       J_dot(i) = j_dot(i) + 4._WP*Omega_orb*q**2/SQRT(R_a*(1+q))*(R_a)**(l+3)* &
            kappa*abs(F)*sin(gamma)*secular_G_4(R_a, e, l, m, k)

    endif

    ! Finish

    return

  end subroutine process_wave

end program gyre_orbit
