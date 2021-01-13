! Module   : gyre_bracket_search
! Purpose  : mode searching (real, bracketing)
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

module gyre_bracket_search

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_discrim_func
  use gyre_ext
  use gyre_mode
  use gyre_num_par
  use gyre_rad_bvp
  use gyre_root
  use gyre_state
  use gyre_status
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  integer, save :: j_m = 0

  ! Access specifiers

  private

  public :: bracket_search

contains

  subroutine bracket_search (bp, omega, omega_min, omega_max, nm_p, process_mode)

    class(r_bvp_t), target, intent(inout) :: bp
    real(WP), intent(in)                  :: omega(:)
    real(WP), intent(in)                  :: omega_min
    real(WP), intent(in)                  :: omega_max
    type(num_par_t), intent(in)           :: nm_p
    interface
       subroutine process_mode (md, n_iter, chi)
         use core_kinds
         use gyre_ext
         use gyre_mode
         type(mode_t), intent(in)  :: md
         integer, intent(in)       :: n_iter
         type(r_ext_t), intent(in) :: chi
       end subroutine process_mode
    end interface

    type(r_state_t)            :: st
    type(r_discrim_func_t)     :: df
    real(WP), allocatable      :: omega_a(:)
    real(WP), allocatable      :: omega_b(:)
    type(r_ext_t), allocatable :: discrim_a(:)
    type(r_ext_t), allocatable :: discrim_b(:)
    integer                    :: c_beg
    integer                    :: c_end
    integer                    :: c_rate
    integer                    :: n_iter
    integer                    :: i
    type(r_ext_t)              :: omega_root
    integer                    :: status
    type(wave_t)               :: wv
    type(mode_t)               :: md
    type(r_ext_t)              :: chi

    ! Set up the discriminant function

    st = r_state_t()
    df = r_discrim_func_t(bp, st, omega_min, omega_max)

    ! Find discriminant root brackets

    call find_brackets_(df, omega, omega_a, omega_b, discrim_a, discrim_b)

    ! Process each bracket to find modes

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'm', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(1X,A3,1X,A4,1X,A7,1X,A6,1X,A6,1X,A15,1X,A15,1X,A10,1X,A6)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    mode_loop : do i = 1, SIZE(omega_a)

       ! Find the discriminant root

       n_iter = 0

       call solve(df, r_ext_t(omega_a(i)), r_ext_t(omega_b(i)), r_ext_t(0._WP), nm_p, &
                  omega_root, status, n_iter=n_iter, n_iter_max=nm_p%n_iter_max, f_rx_a=discrim_a(i), f_rx_b=discrim_b(i))
       if (status /= STATUS_OK) then
          call report_status_(status, 'solve')
          cycle mode_loop
       endif

       ! Construct the mode_t

       j_m = j_m + 1

       st = r_state_t(omega=real(omega_root))

       select type (bp)
       type is (ad_bvp_t)
          wv = wave_t(bp, st, j_m)
       type is (rad_bvp_t)
          wv = wave_t(bp, st, j_m)
       class default
          $ABORT(Invalid bp class)
       end select

       md = mode_t(wv)

       ! Process it

       chi = abs(md%discrim)/max(abs(discrim_a(i)), abs(discrim_b(i)))
       
       call process_mode(md, n_iter, chi)

    end do mode_loop

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(2X,A,1X,F10.3,1X,A/)
    endif

    ! Finish

    return

  contains

    subroutine report_status_ (status, stage_str)

      integer, intent(in)      :: status
      character(*), intent(in) :: stage_str

      ! Report the status

      if (check_log_level('WARN')) then

         write(OUTPUT_UNIT, 100) 'Failed during ', stage_str, ' : ', status_str(status)
100      format(4A)

      endif

      if (check_log_level('INFO')) then

         write(OUTPUT_UNIT, 110) 'n_iter  :', n_iter
110      format(3X,A,1X,I0)

         write(OUTPUT_UNIT, 120) 'omega_a :', omega_a(i)
         write(OUTPUT_UNIT, 120) 'omega_b :', omega_b(i)
120      format(3X,A,1X,E24.16)

      end if

      ! Finish

      return

    end subroutine report_status_
      
  end subroutine bracket_search

  !****

  subroutine find_brackets_ (df, omega, omega_a, omega_b, discrim_a, discrim_b)

    type(r_discrim_func_t), intent(inout)   :: df
    real(WP), intent(in)                    :: omega(:)
    real(WP), allocatable, intent(out)      :: omega_a(:)
    real(WP), allocatable, intent(out)      :: omega_b(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_a(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_b(:)

    integer              :: n_omega
    integer, allocatable :: i_part(:)
    integer              :: c_beg
    integer              :: c_end
    integer              :: c_rate
    integer              :: i
    type(r_ext_t)        :: discrim_i
    integer              :: status
    real(WP)             :: discrim_f(SIZE(omega))
    integer              :: discrim_e(SIZE(omega))
    $if ($MPI)
    integer              :: p
    $endif
    type(r_ext_t)        :: discrim(SIZE(omega))
    integer              :: n_brack
    integer              :: i_brack(SIZE(omega))

    ! Calculate the discriminant on the omega abscissa

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Evaluating discriminant'
100    format(A)
    endif

    n_omega = SIZE(omega)

    allocate(i_part(MPI_SIZE+1))

    call partition_tasks(n_omega, 1, i_part)

    call SYSTEM_CLOCK(c_beg, c_rate)

    discrim_loop: do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       call df%eval(r_ext_t(omega(i)), discrim_i, status)

       discrim_f(i) = FRACTION(discrim_i)
       discrim_e(i) = EXPONENT(discrim_i)

       if (check_log_level('DEBUG')) then
          write(OUTPUT_UNIT, 110) omega(i), discrim_f(i), discrim_e(i)
110       format(2X,E24.16,2X,F19.16,2X,I7)
       endif

    end do discrim_loop

    $if ($MPI)

    do p = 1,MPI_SIZE
       call bcast_seq(discrim_f, i_part(p), i_part(p+1)-1, p-1)
       call bcast_seq(discrim_e, i_part(p), i_part(p+1)-1, p-1)
    end do

    $endif

    discrim = scale(r_ext_t(discrim_f), discrim_e)

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
120    format(2X,A,F10.3,1X,A/)
    endif

    ! Find root brackets

    n_brack = 0

    bracket_loop : do i = 1, n_omega-1

       if (discrim(i)*discrim(i+1) <= r_ext_t(0._WP)) then
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

  end subroutine find_brackets_

end module gyre_bracket_search
