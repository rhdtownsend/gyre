! Program  : gyre_response
! Purpose  : tidal response code
!
! Copyright 2019 Rich Townsend & The GYRE Team
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

program gyre_response

  ! Uses

  use core_kinds, only : WP
  use core_hgroup
  use core_parallel
  use core_system

  use gyre_context
  use gyre_constants
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
  complex(WP), allocatable       :: xi_r(:,:,:)
  complex(WP), allocatable       :: lag_L(:,:,:)
  integer                        :: i
  integer                        :: l
  integer                        :: m
  integer                        :: k
  character(LEN=FILENAME_LEN)    :: response_file
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

  ! Allocate displacement and luminosity perturbations arrays

  allocate(xi_r(0:td_p(1)%l_max,-td_p(1)%l_max:td_p(1)%l_max,0:td_p(1)%k_max))
  allocate(lag_L(0:td_p(1)%l_max,-td_p(1)%l_max:td_p(1)%l_max,0:td_p(1)%k_max))

  ! Loop over orbital frequencies

  n_Omega_orb = SIZE(Omega_orb)

  Omega_orb_loop : do i = 1, n_Omega_orb

     ! Initialize the perturbation arrays

     xi_r = 0._WP
     lag_L = 0._WP

     ! Add in contributions from each tidal component

     call eval_tide(ml, process_wave, Omega_orb(i), os_p(1), nm_p(1), gr_p(1), td_p(1))

     ! Perform sanity checks

     if (.FALSE.) then

        do l = 2, td_p(1)%l_max
           do m = -l, l
              do k = -td_p(1)%k_max, td_p(1)%k_max
                 if (ABS(xi_r(l,m,k)) /= 0._WP) then
                    print 125, l, m, k, ABS(xi_r(l,m,k) - (-1)**m*CONJG(xi_r(l,-m,-k)))/ABS(xi_r(l,m,k))
125                 format(I3,1X,I3,1X,I3,F14.10)
                 endif
              end do
           end do
        end do

     end if
        
     ! Write out results

     write(response_file, 130) 'response.', i, '.h5'
130  format(A,I3.3,A)

     hg = hgroup_t(response_file, CREATE_FILE)

     call write_attr(hg, 'Omega_orb', Omega_orb(i))
     call write_attr(hg, 'Omega_rot', ml_p%Omega_rot)

     call write_attr(hg, 'l_max', td_p(1)%l_max)
     call write_attr(hg, 'k_max', td_p(1)%k_max)

     call write_attr(hg, 'e', td_p(1)%e)

     call write_dset(hg, 'xi_r', xi_r)
     call write_dset(hg, 'lag_L', lag_L)

     call hg%final()

  end do Omega_orb_loop

  ! Clean up

  deallocate(ml)

  ! Finish

  call final_parallel()

contains

  subroutine process_wave (wv, k)

    type(wave_t), intent(in) :: wv
    integer, intent(in)      :: k

    integer :: l
    integer :: m

    l = wv%md_p%l
    m = wv%md_p%m

    ! Store the displacement

    xi_r(l,m,k) = wv%xi_r(wv%k_ref)

    ! Store the luminosity perturbation

    lag_L(l,m,k) = wv%lag_L(wv%k_ref)

    ! Finish

    return

  end subroutine process_wave

end program gyre_response
