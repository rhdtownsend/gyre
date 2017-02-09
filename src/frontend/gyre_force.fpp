! Program  : gyre_force
! Purpose  : forced oscillation code
!
! Copyright 2016-2017 Rich Townsend
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
  use gyre_freq
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode
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
  use gyre_util
  use gyre_version

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
  real(WP)                      :: q
  real(WP)                      :: c_lmk
  integer                       :: n_omega
  real(WP), allocatable         :: omega(:)
  integer                       :: n_P
  real(WP), allocatable         :: P(:)
  real(WP)                      :: M_pri
  real(WP)                      :: M_sec
  real(WP)                      :: R_pri
  integer                       :: i
  type(osc_par_t)               :: os_p_sel
  type(num_par_t)               :: nm_p_sel
  type(grid_par_t)              :: gr_p_sel
  type(grid_t)                  :: gr
  integer                       :: j
  class(r_bvp_t), allocatable   :: bp_ad
  class(c_bvp_t), allocatable   :: bp_nad

  namelist /force/ q, c_lmk, n_omega, omega, n_P, P

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

  ! Initialize binary masses

  select type (ml)
  class is (evol_model_t)
     M_pri = ml%M_star
     M_sec = q*M_pri
     R_pri = ml%R_star
  class default
     $ABORT(Invalid class)
  end select

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

     ! Create the scaffold grid (used in setting up the frequency array)

     gr = grid_t(ml%grid(), gr_p_sel%x_i, gr_p_sel%x_o)

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

     ! Create the full grid

     gr = grid_t(ml, omega, gr_p_sel, md_p(i), os_p_sel)

     ! Find modes

     if (os_p_sel%nonadiabatic) then

        allocate(bp_nad, SOURCE=nad_bvp_t(ml, gr, md_p(i), nm_p_sel, os_p_sel))

        call scan_force_c(bp_nad, omega, P)

        deallocate(bp_nad)

     else

        if (md_p(i)%l == 0 .AND. os_p_sel%reduce_order) then
           allocate(bp_ad, SOURCE=rad_bvp_t(ml, gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bvp_t(ml, gr, md_p(i), nm_p_sel, os_p_sel))
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

    type(${T}_bvp_t), intent(inout) :: bp
    real(WP), intent(in)            :: omega(:)
    real(WP), intent(in)            :: P(:)

    integer   :: n_omega
    integer   :: j
    real(WP)  :: a
    real(WP)  :: eps_T
    real(WP)  :: alpha_fc
    $TYPE(WP) :: w_i(bp%n_i)
    $TYPE(WP) :: w_o(bp%n_i)
    $TYPE(WP) :: y(bp%n_e,bp%n_k)

    character(64) :: filename
    integer       :: res_unit
    integer       :: sol_unit
    integer       :: k

    $CHECK_BOUNDS(SIZE(P),SIZE(omega))

    ! Scan over frequencies

    open(NEWUNIT=res_unit, FILE='response.dat', STATUS='replace')

    n_omega = SIZE(omega)

    omega_loop : do j = 1, n_omega

       ! Set up binary parameters

       a = (G_GRAVITY*(M_pri + M_sec)*P(j)**2/(4.*PI**2))**(1._WP/3._WP)

       eps_T = (R_pri/a)**3*(M_sec/M_pri)

       alpha_fc = eps_T*(2*md_p(i)%l+1)*c_lmk

       ! Set up the inhomogeneous boundary terms

       w_i = 0._WP

       w_o = 0._WP
       w_o(2) = -alpha_fc

       ! Build and solve the linear system

       $if($T eq 'c')
       call bp%build(CMPLX(omega(j), KIND=WP))
       $else
       call bp%build(omega(j))
       $endif

       y = bp%soln_vec_inhom(w_i, w_o)

       ! Write out the tidal response

       write(res_unit, 100) omega(j), P(j), y(1,bp%n_k)
100    format(999E16.8)

       ! Write out the solution data

       write(filename, 110) 'forced.', j, '.txt'
110    format(A,I3.3,A)

       open(NEWUNIT=sol_unit, FILE=filename, STATUS='REPLACE')

       do k = 1,bp%n_k
          write(sol_unit, 120) gr%pt(k)%x, y(:,k)
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
