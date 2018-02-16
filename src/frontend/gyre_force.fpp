! Program  : gyre_force
! Purpose  : forced oscillation code
!
! Copyright 2016-2018 Rich Townsend
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

program gyre_force

  ! Uses

  use core_kinds, only : WP
  use core_parallel
  use core_system

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_ext
  use gyre_evol_model
  use gyre_context
  use gyre_freq
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
  use gyre_rad_bvp
  use gyre_scan_par
  use gyre_search
  use gyre_state
  use gyre_util
  use gyre_version
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: D_0 = 1024

  ! Variables

  character(:), allocatable     :: filename
  integer                       :: unit
  type(model_par_t)             :: ml_p
  type(mode_par_t), allocatable :: md_p(:)
  type(osc_par_t), allocatable  :: os_p(:)
  type(num_par_t), allocatable  :: nm_p(:)
  type(grid_par_t), allocatable :: gr_p(:)
  class(model_t), pointer       :: ml => null()
  type(context_t), pointer      :: cx(:) => null()
  real(WP)                      :: F
  integer                       :: n_omega
  real(WP), allocatable         :: omega(:)
  integer                       :: n_P
  real(WP), allocatable         :: P(:)
  integer                       :: i
  type(osc_par_t)               :: os_p_sel
  type(num_par_t)               :: nm_p_sel
  type(grid_par_t)              :: gr_p_sel
  type(grid_t)                  :: gr
  integer                       :: j
  class(r_bvp_t), allocatable   :: bp_ad
  class(c_bvp_t), allocatable   :: bp_nad

  namelist /force/ F, n_omega, omega, n_P, P

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_force <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_force ['//VERSION//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110  format(A,1X,A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

  endif

  ! Read the namelist file

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_mode_par(unit, md_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)

  allocate(omega(D_0))
  allocate(P(D_0))

  n_omega = 0
  n_P = 0

  rewind(unit)
  read(unit, NML=force)

  $ASSERT(n_omega > 0 .NEQV. n_P > 0,One or other of n_omega/n_P must be non-zero)

  $ASSERT(n_omega <= D_0,Frequency array too short)
  $ASSERT(n_P <= D_0,Period array too short)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '=')
  endif

  ml => model_t(ml_p)

  ! Allocate the contexts array (will be initialized later on)

  allocate(cx(SIZE(md_p)))

  ! Loop through md_p

  md_p_loop : do i = 1, SIZE(md_p)

     if (check_log_level('INFO')) then

        write(OUTPUT_UNIT, 100) form_header('Mode Search', '=')

        write(OUTPUT_UNIT, 100) 'Mode parameters'

        write(OUTPUT_UNIT, 130) 'l :', md_p(i)%l
        write(OUTPUT_UNIT, 130) 'm :', md_p(i)%m
130     format(3X,A,1X,I0)

        write(OUTPUT_UNIT, *)

     endif

     ! Select parameters according to tags

     call select_par(os_p, md_p(i)%tag, os_p_sel)
     call select_par(nm_p, md_p(i)%tag, nm_p_sel)
     call select_par(gr_p, md_p(i)%tag, gr_p_sel)

     ! Set up the context

     cx(i) = context_t(ml, gr_p_sel, md_p(i), os_p_sel)

     ! Set up the frequency/period arrays

     if (n_omega == 0) then
        do j = 1, n_P
           omega(j) = omega_from_freq(1._WP/P(j), ml, gr%pt(1), gr%pt(gr%n_k), 'HZ', 'INERTIAL', md_p(i), os_p_sel)
        end do
        n_omega = n_P
     elseif (n_P == 0) then
        do j = 1, n_omega
           P(j) = 1._WP/freq_from_omega(omega(j), ml, gr%pt(1), gr%pt(gr%n_k), 'HZ', 'INERTIAL', md_p(i), os_p_sel)
        end do
        n_P = n_omega
     endif

     omega = omega(:n_omega)
     P = P(:n_P)

     ! Create the grid

     gr = grid_t(cx(i), omega, gr_p_sel)

     ! Find modes

     if (os_p_sel%nonadiabatic) then

        allocate(bp_nad, SOURCE=nad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))

        call scan_force_c(bp_nad, omega, P)

        deallocate(bp_nad)

     else

        if (md_p(i)%l == 0 .AND. os_p_sel%reduce_order) then
           allocate(bp_ad, SOURCE=rad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))
        endif

        call scan_force_r(bp_ad, omega, P)

        deallocate(bp_ad)

     end if

  end do md_p_loop

  ! Clean up

  deallocate(ml)

  ! Finish

  close(unit)

  call final_parallel()

contains

  $define $SCAN_FORCE $sub

  $local $T $1
  $local $TYPE $2

  subroutine scan_force_${T} (bp, omega, P)

    class(${T}_bvp_t), intent(inout) :: bp
    real(WP), intent(in)             :: omega(:)
    real(WP), intent(in)             :: P(:)

    integer            :: n_omega
    integer            :: j
    $TYPE(WP)          :: w_i(bp%n_i)
    $TYPE(WP)          :: w_o(bp%n_o)
    type(${T}_state_t) :: st

    character(64) :: filename
    type(wave_t)  :: wv
    integer       :: res_unit
    integer       :: sol_unit
    integer       :: i_
    integer       :: k_
    real(WP)      :: tau

    $CHECK_BOUNDS(SIZE(P),SIZE(omega))

    ! Scan over frequencies

    open(NEWUNIT=res_unit, FILE='response.dat', STATUS='replace')

    n_omega = SIZE(omega)

    omega_loop : do j = 1, n_omega

       ! Set up the inhomogeneous boundary terms

       w_i = 0._WP
         
       w_o = 0._WP
       w_o(2) = F
         
       ! Solve for the wave function

       $if($T eq 'c')
       st = c_state_t(CMPLX(omega(j), KIND=WP), 0._WP)
       select type (bp)
       type is (nad_bvp_t)
          wv = wave_t(bp, st, w_i, w_o)
       class default
          $ABORT(Invalid bp class)
       end select
       $else
       st = r_state_t(omega(j))
       select type (bp)
       type is (ad_bvp_t)
          wv = wave_t(bp, st, w_i, w_o)
       type is (rad_bvp_t)
          wv = wave_t(bp, st, w_i, w_o)
       class default
          $ABORT(Invalid bp class)
       end select
       $endif

       ! Calculate the torque expected

       tau = -0.5_WP*md_p(i)%m*AIMAG(F*CONJG(wv%y_i(3,bp%n_k)))

       ! Write out the response

       write(res_unit, 100) omega(j), P(j), wv%y_i(1,bp%n_k), wv%tau_ss(), tau
100    format(999E16.8)

       ! Write out the solution data

       write(filename, 110) 'forced.', j, '.txt'
110    format(A,I3.3,A)

       open(NEWUNIT=sol_unit, FILE=filename, STATUS='REPLACE')

       do k_ = 1,bp%n_k
          write(sol_unit, 120) wv%gr%pt(k_)%x, (wv%y_i(i_,k_),i_=1,6), wv%dtau_dx_ss(k_)
120       format(999(1X,E16.8))
       enddo

       close(sol_unit)

    end do omega_loop

    close(res_unit)

    ! Finish

    return

  end subroutine scan_force_${T}

  $endsub

  $SCAN_FORCE(r,real)
  $SCAN_FORCE(c,complex)

end
