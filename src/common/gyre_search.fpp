! Module   : gyre_search
! Purpose  : mode searching
!
! Copyright 2013 Rich Townsend
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

module gyre_search

  ! Uses

  use core_kinds
  use gyre_constants
  use core_order
  use core_parallel

  use gyre_bvp
  use gyre_model
  use gyre_modepar
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_mode
  use gyre_grid
  use gyre_ext_arith
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_scan
  public :: scan_search
  public :: prox_search

contains

  subroutine build_scan (sp, ml, mp, op, gp, x_in, omega)

    type(scanpar_t), intent(in)        :: sp(:)
    class(model_t), intent(in)         :: ml
    type(modepar_t), intent(in)        :: mp
    type(oscpar_t), intent(in)         :: op
    type(gridpar_t), intent(in)        :: gp(:)
    real(WP), allocatable, intent(in)  :: x_in(:)
    real(WP), allocatable, intent(out) :: omega(:)

    real(WP)              :: x_i
    real(WP)              :: x_o
    integer               :: i
    real(WP)              :: omega_min
    real(WP)              :: omega_max
    integer               :: j
    integer               :: n_omega
    real(WP)              :: omega_cutoff_lo
    real(WP)              :: omega_cutoff_hi

    $ASSERT(SIZE(sp) >=1,Empty scanpars)

    ! Determine the grid range

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Building omega grid'
100    format(A)
    endif

    call grid_range(gp, ml, mp, x_in, x_i, x_o)

    ! Loop through scanpars

    allocate(omega(0))

    sp_loop : do i = 1,SIZE(sp)

       ! Determine the frequency range

       omega_min = sp(i)%freq_min/freq_scale(ml, mp, op, x_o, sp(i)%freq_units)
       omega_max = sp(i)%freq_max/freq_scale(ml, mp, op, x_o, sp(i)%freq_units)

       ! Add points to the frequency grid

       select case(sp(i)%grid_type)
       case('LINEAR')
          omega = [omega,(((sp(i)%n_freq-j)*omega_min + (j-1)*omega_max)/(sp(i)%n_freq-1), j=1,sp(i)%n_freq)]
       case('INVERSE')
          omega = [omega,((sp(i)%n_freq-1)/((sp(i)%n_freq-j)/omega_min + (j-1)/omega_max), j=1,sp(i)%n_freq)]
       case default
          $ABORT(Invalid grid_type)
       end select

    end do sp_loop

    n_omega = SIZE(omega)

    $ASSERT(n_omega > 2,At least two frequency points required)

    ! Sort the frequencies

    omega = omega(sort_indices(omega))

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 110) 'omega points :', n_omega
110    format(2X,A,1X,I0)

       write(OUTPUT_UNIT, 120) 'omega range  :', MINVAL(omega), '->',  MAXVAL(omega)
120    format(2X,A,1X,E24.16,1X,A,1X,E24.16)

    endif

    ! Perform checks

    if (check_log_level('WARN')) then

       call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)

       if (omega(1) < omega_cutoff_lo) then
          write(OUTPUT_UNIT, 100) '!!! WARNING: omega extends below atmospheric gravity cutoff frequency'
130       format(2X,A)
       end if

       if (omega(n_omega) > omega_cutoff_hi) then
          write(OUTPUT_UNIT, 100) '!!! WARNING: omega extends above atmospheric acoustic cutoff frequency'
       end if

    endif

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine build_scan

!****

  subroutine scan_search (bp, omega, process_mode)

    class(bvp_t), intent(inout) :: bp
    real(WP), intent(in)        :: omega(:)
    interface
       subroutine process_mode (md)
         use gyre_mode
         type(mode_t), intent(in) :: md
       end subroutine process_mode
    end interface

    real(WP), allocatable         :: omega_a(:)
    real(WP), allocatable         :: omega_b(:)
    type(ext_real_t), allocatable :: discrim_a(:)
    type(ext_real_t), allocatable :: discrim_b(:)
    integer                       :: n_brack
    integer                       :: c_beg
    integer                       :: c_end
    integer                       :: c_rate
    type(mode_t)                  :: md
    integer                       :: i

    ! Find discriminant root brackets

    call find_brackets(bp, omega, omega_a, omega_b, discrim_a, discrim_b)

    ! Process each bracket to find modes

    if(check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter', 'n'
110    format(4(2X,A8),3(2X,A24),2X,A6,2X,A7)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    mode_loop : do i = 1, SIZE(omega_a)

       ! Find the mode

       md = bp%mode(CMPLX([omega_a(i),omega_b(i)], KIND=WP), &
                    ext_complex_t([discrim_a(i),discrim_b(i)]), use_real=.TRUE.)
 
       if (md%n_pg < md%mp%n_pg_min .OR. md%n_pg > md%mp%n_pg_max) cycle mode_loop

       ! Process it

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 120) md%mp%l, md%n_pg, md%n_p, md%n_g, md%omega, real(md%chi), md%n_iter, md%n
120       format(4(2X,I8),3(2X,E24.16),2X,I6,2X,I7)
       endif

       call process_mode(md)

       $if($GFORTRAN_PR57922)
       call md%final()
       $endif

    end do mode_loop

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(2X,A,1X,F10.3,1X,A)
    endif

    ! Finish

    return

  end subroutine scan_search

!****

  subroutine prox_search (bp, md_in, process_mode)

    class(bvp_t), target, intent(inout) :: bp
    type(mode_t), intent(in)            :: md_in(:)
    interface
       subroutine process_mode (md)
         use gyre_mode
         type(mode_t), intent(in) :: md
       end subroutine process_mode
    end interface
    
    integer      :: c_beg
    integer      :: c_end
    integer      :: c_rate
    integer      :: i
    complex(WP)  :: omega_a
    complex(WP)  :: omega_b
    type(mode_t) :: md

    ! Process each initial mode to find a proximate mode

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter', 'n'
110    format(4(2X,A8),3(2X,A24),2X,A6,2X,A7)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    mode_loop : do i = 1, SIZE(md_in)

       ! Set up initial guesses

       omega_a = md_in(i)%omega*CMPLX(1._WP,  SQRT(EPSILON(0._WP)), WP)
       omega_b = md_in(i)%omega*CMPLX(1._WP, -SQRT(EPSILON(0._WP)), WP)

       ! Find the mode

       md = bp%mode([omega_a,omega_b])

       if (md%n_pg < md%mp%n_pg_min .OR. md%n_pg > md%mp%n_pg_max) cycle mode_loop

       ! Process it

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 120) md%mp%l, md%n_pg, md%n_p, md%n_g, md%omega, real(md%chi), md%n_iter, md%n
120       format(4(2X,I8),3(2X,E24.16),2X,I6,2X,I7)
       endif

       call process_mode(md)

       $if($GFORTRAN_PR57922)
       call md%final()
       $endif

    end do mode_loop

    call SYSTEM_CLOCK(c_end)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(2X,A,1X,F10.3,1X,A)
    endif

    ! Finish

    return

  end subroutine prox_search

!****

  subroutine find_brackets (bp, omega, omega_a, omega_b, discrim_a, discrim_b)

    class(bvp_t), target, intent(inout)        :: bp
    real(WP), intent(in)                       :: omega(:)
    real(WP), allocatable, intent(out)         :: omega_a(:)
    real(WP), allocatable, intent(out)         :: omega_b(:)
    type(ext_real_t), allocatable, intent(out) :: discrim_a(:)
    type(ext_real_t), allocatable, intent(out) :: discrim_b(:)

    integer          :: n_omega
    integer          :: c_beg
    integer          :: c_end
    integer          :: c_rate
    integer          :: i
    type(ext_real_t) :: discrim(SIZE(omega))
    integer          :: n_brack
    integer          :: i_brack(SIZE(omega))

    ! Calculate the discriminant on the omega abscissa

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Root bracketing'
100    format(A)
    endif

    n_omega = SIZE(omega)

    call SYSTEM_CLOCK(c_beg, c_rate)

    discrim_loop : do i = 1, n_omega

       discrim(i) = ext_real_t(bp%discrim(CMPLX(omega(i), KIND=WP), use_real=.TRUE.))

       if(check_log_level('DEBUG')) then
          write(OUTPUT_UNIT, 110) omega(i), fraction(discrim(i)), exponent(discrim(i))
110       format(2X,E24.16,2X,F19.16,2X,I7)
       endif

    end do discrim_loop

    call SYSTEM_CLOCK(c_end)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
120    format(2X,A,F10.3,1X,A/)
    endif

    ! Find root brackets

    n_brack = 0

    bracket_loop : do i = 1, n_omega-1

       if(discrim(i)*discrim(i+1) <= ext_real_t(0._WP)) then
          n_brack = n_brack + 1
          i_brack(n_brack) = i
       end if

    end do bracket_loop

    ! Set up the bracket frequencies

    omega_a = omega(i_brack(:n_brack))
    omega_b = omega(i_brack(:n_brack)+1)

    discrim_a = discrim(i_brack(:n_brack))
    discrim_b = discrim(i_brack(:n_brack)+1)

    ! Finish

    return

  end subroutine find_brackets

end module gyre_search
