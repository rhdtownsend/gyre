! Module   : gyre_minmod_search
! Purpose  : mode searching (complex, minmod)
!
! Copyright 2013-2021 Rich Townsend & The GYRE Team
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

module gyre_minmod_search

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_bvp
  use gyre_discrim
  use gyre_ext
  use gyre_num_par
  use gyre_min
  use gyre_prox_search
  use gyre_state
  use gyre_status
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  integer, save :: j_m = 0

  ! Access specifiers

  private

  public :: minmod_search

contains

  subroutine minmod_search (bp, omega, omega_min, omega_max, nm_p, process_mode)

    class(c_bvp_t), intent(inout) :: bp
    real(WP), intent(in)          :: omega(:)
    real(WP), intent(in)          :: omega_min
    real(WP), intent(in)          :: omega_max
    type(num_par_t), intent(in)   :: nm_p
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

    real(WP), allocatable      :: omega_a(:)
    real(WP), allocatable      :: omega_b(:)
    real(WP), allocatable      :: omega_c(:)
    type(r_ext_t), allocatable :: discrim_a(:)
    type(r_ext_t), allocatable :: discrim_b(:)
    type(r_ext_t), allocatable :: discrim_c(:)
    integer                    :: n_in
    complex(WP), allocatable   :: omega_in(:)
    integer, allocatable       :: j_in(:)
    type(r_ext_t)              :: x_a
    type(r_ext_t)              :: x_b
    type(r_ext_t)              :: x_c
    type(r_ext_t)              :: f_x_a
    type(r_ext_t)              :: f_x_b
    type(r_ext_t)              :: f_x_c
    integer                    :: i
    integer                    :: n_iter
    type(r_ext_t)              :: x_min
    integer                    :: status

    ! Find discriminant minmod brackets

    call find_brackets_(bp, omega, omega_a, omega_b, omega_c, discrim_a, discrim_b, discrim_c)

    ! Convert brackets into minmod frequencies

    n_in = SIZE(omega_a)

    allocate(omega_in(n_in))
    allocate(j_in(n_in))

    in_loop : do i = 1, n_in

       ! Solve for the discriminant minimum

       n_iter = 0

       x_a = r_ext_t(omega_a(i))
       x_b = r_ext_t(omega_b(i))
       x_c = r_ext_t(omega_c(i))
       f_x_a = discrim_a(i)
       f_x_b = discrim_b(i)
       f_x_c = discrim_c(i)

       call solve(eval_discrim_x_, x_a, x_b, x_c, r_ext_t(0._WP), nm_p, &
                  x_min, status, n_iter=n_iter, n_iter_max=nm_p%n_iter_max, &
                  f_x_a=f_x_a, f_x_b=f_x_b, f_x_c=f_x_c)
       if (status /= STATUS_OK) then
          call report_status_(status, 'solve')
          cycle in_loop
       endif

       j_m = j_m + 1

       omega_in(i) = real(x_min)
       j_in(i) = j_m

    end do in_loop

    ! Search for modes

    call prox_search(bp, omega_in, j_in, omega_min, omega_max, nm_p, process_mode)

    ! Finish

    return

  contains

    subroutine eval_discrim_x_ (x, discrim, status)

      type(r_ext_t), intent(in)  :: x
      type(r_ext_t), intent(out) :: discrim
      integer, intent(out)       :: status

      real(WP)      :: omega_r
      complex(WP)   :: omega
      type(c_ext_t) :: discrim_c

      ! Evaluate the frequency from x
      
      omega_r = real(x)
      omega = CMPLX(omega_r, KIND=WP)

      ! Evaluate the discriminant

      if (omega_r >= omega_min .AND. omega_r <= omega_max) then

         call eval_discrim(bp, c_state_t(omega, omega_r=omega_r), discrim_c)

         discrim = abs(discrim_c)

         status = STATUS_OK

      else

         status = STATUS_OMEGA_DOMAIN

      endif
      
      ! Finish

      return

    end subroutine eval_discrim_x_

    !****

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
         write(OUTPUT_UNIT, 120) 'omega_c :', omega_c(i)
120      format(3X,A,1X,E24.16)

      end if

      ! Finish

      return

    end subroutine report_status_
      
  end subroutine minmod_search

  !****

  subroutine find_brackets_ (bp, omega, omega_a, omega_b, omega_c, discrim_a, discrim_b, discrim_c)

    class(c_bvp_t), intent(inout)           :: bp
    real(WP), intent(in)                    :: omega(:)
    real(WP), allocatable, intent(out)      :: omega_a(:)
    real(WP), allocatable, intent(out)      :: omega_b(:)
    real(WP), allocatable, intent(out)      :: omega_c(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_a(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_b(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_c(:)

    integer              :: n_omega
    integer, allocatable :: i_part(:)
    integer              :: c_beg
    integer              :: c_end
    integer              :: c_rate
    integer              :: i
    type(c_ext_t)        :: discrim_i
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

       call eval_discrim(bp, c_state_t(CMPLX(omega(i), KIND=WP), omega_r=omega(i)), discrim_i)

       discrim_f(i) = FRACTION(abs(discrim_i))
       discrim_e(i) = EXPONENT(abs(discrim_i))

       if (check_log_level('DEBUG')) then
          write(OUTPUT_UNIT, 110) omega(i), fraction(discrim(i)), exponent(discrim(i))
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

    ! Find minimum brackets

    n_brack = 0

    bracket_loop : do i = 2, n_omega-1

       if (discrim(i) < discrim(i-1) .AND. discrim(i) < discrim(i+1)) then
          n_brack = n_brack + 1
          i_brack(n_brack) = i
       end if

    end do bracket_loop

    ! Set up the bracket frequencies

    omega_a = omega(i_brack(:n_brack)-1)
    omega_b = omega(i_brack(:n_brack))
    omega_c = omega(i_brack(:n_brack)+1)

    discrim_a = discrim(i_brack(:n_brack)-1)
    discrim_b = discrim(i_brack(:n_brack))
    discrim_c = discrim(i_brack(:n_brack)+1)

    ! Finish

    return

  end subroutine find_brackets_

end module gyre_minmod_search
