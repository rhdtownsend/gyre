! Program  : gyre_force
! Purpose  : forced oscillation code
!
! Copyright 2016 Rich Townsend
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

  use gyre_ad_bep
  use gyre_bep
  use gyre_constants
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bep
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
  use gyre_rad_bep
  use gyre_scan_par
  use gyre_search
  use gyre_util
  use gyre_version

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable     :: filename
  integer                       :: unit
  type(model_par_t)             :: ml_p
  type(mode_par_t), allocatable :: md_p(:)
  type(osc_par_t), allocatable  :: os_p(:)
  type(num_par_t), allocatable  :: nm_p(:)
  type(grid_par_t), allocatable :: gr_p(:)
  type(scan_par_t), allocatable :: sc_p(:)
  class(model_t), pointer       :: ml => null()
  integer                       :: i
  type(osc_par_t)               :: os_p_sel
  type(num_par_t)               :: nm_p_sel
  type(grid_par_t)              :: gr_p_sel
  type(scan_par_t), allocatable :: sc_p_sel(:)
  type(grid_t)                  :: gr
  real(WP), allocatable         :: omega(:)
  integer                       :: s
  integer                       :: k_i
  integer                       :: k_o
  class(r_bep_t), allocatable   :: bp_ad
  class(c_bep_t), allocatable   :: bp_nad

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
  call read_scan_par(unit, sc_p)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '=')
  endif

  ml => model_t(ml_p)

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
     call select_par(sc_p, md_p(i)%tag, sc_p_sel)

     ! Create the scaffold grid (used in setting up the frequency array)

     gr = grid_t(ml%grid(), gr_p_sel%x_i, gr_p_sel%x_o)

     ! Set up the frequency array

     call build_scan(ml, gr, md_p(i), os_p_sel, sc_p_sel, omega)

     call check_scan(ml, gr, omega, md_p(i), os_p_sel)

     ! Create the full grid

     if (check_log_level('INFO')) then
        write(OUTPUT_UNIT, 100) 'Building x grid'
     endif

     gr = grid_t(ml, omega, gr_p_sel, md_p(i), os_p_sel)

     if (check_log_level('INFO')) then

        seg_loop : do s = gr%s_i(), gr%s_o()
           k_i = gr%k_i(s)
           k_o = gr%k_o(s)
           write(OUTPUT_UNIT, 140) 'segment', s, ':', k_o-k_i+1, 'points, x range', gr%pt(k_i)%x, '->', gr%pt(k_o)%x
140        format(3X,A,1X,I0,1X,A,1X,I0,1X,A,1X,F6.4,1X,A,1X,F6.4)
        end do seg_loop

        write(OUTPUT_UNIT, *)
        
     end if

     ! Find modes

     if (os_p_sel%nonadiabatic) then

        $ASSERT(ml%nonad_cap(),Model does not have capability for nonadibatic calculations)

        allocate(bp_nad, SOURCE=nad_bep_t(ml, gr, omega, md_p(i), nm_p_sel, os_p_sel))

        call scan_force_c(bp_nad, omega)

        deallocate(bp_nad)

     else

        if (md_p(i)%l == 0 .AND. os_p_sel%reduce_order) then
           allocate(bp_ad, SOURCE=rad_bep_t(ml, gr, omega, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bep_t(ml, gr, omega, md_p(i), nm_p_sel, os_p_sel))
        endif

        call scan_force_r(bp_ad, omega)

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

  subroutine scan_force_${T} (bp, omega)

    type(${T}_bep_t), intent(inout) :: bp
    real(WP), intent(in)            :: omega(:)

    $TYPE(WP) :: w_nz(bp%n_e)
    integer   :: n_omega
    integer   :: i
    $TYPE(WP) :: y(bp%n_e,bp%n_k)

    ! Set up the inhomogeneous boundary term

    w_nz = 0._WP

    $if ($T eq 'c')
    w_nz(5) = 1._WP
    $else
    w_nz(4) = 1._WP
    $endif

    ! Scan over frequencies

    n_omega = SIZE(omega)

    omega_loop : do i = 1, n_omega

       ! Solve the linear system

       $if($T eq 'c')
       call bp%solve_forced(CMPLX(omega(i), KIND=WP), w_nz, y)
       $else
       call bp%solve_forced(omega(i), w_nz, y)
       $endif

       ! Print out the tidal response

       print 100, omega(i), y(1,bp%n_k)
100    format(999E16.8)

    end do omega_loop

    ! Finish

    return

  end subroutine scan_force_${T}

  $endsub

  $SCAN_FORCE(r,real)
  $SCAN_FORCE(c,complex)

end
