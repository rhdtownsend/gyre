! Module   : gyre_prox_search
! Purpose  : mode searching (complex, proximity)
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

module gyre_prox_search

  ! Uses

  use core_kinds
  use gyre_constants

  use gyre_bvp
  use gyre_discrim_func
  use gyre_ext
  use gyre_math
  use gyre_mode
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_root
  use gyre_state
  use gyre_status
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface prox_search
     module procedure prox_search_1_
     module procedure prox_search_2_
  end interface prox_search

  ! Access specifiers

  private

  public :: prox_search

contains

  subroutine prox_search_1_ (bp, omega_in, j_in, omega_min, omega_max, nm_p, process_mode)

    class(c_bvp_t), target, intent(inout) :: bp
    complex(WP), intent(in)               :: omega_in(:)
    integer, intent(in)                   :: j_in(:)
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

    integer                  :: n_in
    complex(WP), allocatable :: omega_in_a(:)
    complex(WP), allocatable :: omega_in_b(:)
    complex(WP)              :: domega
    integer                  :: i

    $CHECK_BOUNDS(SIZE(j_in),SIZE(omega_in))

    ! Convert initial frequencies into pairs

    n_in = SIZE(omega_in)

    allocate(omega_in_a(n_in))
    allocate(omega_in_b(n_in))

    in_loop : do i = 1, n_in

       domega = sqrt(EPSILON(0._WP))*abs(omega_in(i))*CMPLX(0._WP, 1._WP, KIND=WP)

       omega_in_a(i) = omega_in(i) + domega
       omega_in_b(i) = omega_in(i) - domega

    end do in_loop

    ! Search for modes
    
    call prox_search_2_(bp, omega_in_a, omega_in_b, j_in, omega_min, omega_max, nm_p, process_mode)

    ! Finish

    return

  end subroutine prox_search_1_

  !****

  subroutine prox_search_2_ (bp, omega_in_a, omega_in_b, j_in, omega_min, omega_max, nm_p, process_mode)

    class(c_bvp_t), target, intent(inout) :: bp
    complex(WP), intent(in)               :: omega_in_a(:)
    complex(WP), intent(in)               :: omega_in_b(:)
    integer, intent(in)                   :: j_in(:)
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

    complex(WP), allocatable :: omega_def(:)
    integer                  :: c_beg
    integer                  :: c_end
    integer                  :: c_rate
    integer                  :: i
    type(c_ext_t)            :: omega_a
    type(c_ext_t)            :: omega_b
    real(WP)                 :: omega_r
    integer                  :: n_iter
    integer                  :: n_iter_def
    type(c_state_t)          :: st
    type(c_discrim_func_t)   :: df
    type(c_state_t)          :: st_def
    type(c_discrim_func_t)   :: df_def
    integer                  :: status
    type(c_ext_t)            :: discrim_a
    type(c_ext_t)            :: discrim_b
    type(c_ext_t)            :: discrim_a_rev
    type(c_ext_t)            :: discrim_b_rev
    type(c_ext_t)            :: omega_root
    type(wave_t)             :: wv
    type(mode_t)             :: md
    type(r_ext_t)            :: chi

    $CHECK_BOUNDS(SIZE(j_in),SIZE(omega_in_a))
    $CHECK_BOUNDS(SIZE(omega_in_b),SIZE(omega_in_a))

    ! Initialize the frequency deflation array

    allocate(omega_def(0))

    ! Process each input frequency pair to find modes

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'm', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(1X,A3,1X,A4,1X,A7,1X,A6,1X,A6,1X,A15,1X,A15,1X,A10,1X,A6)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    in_loop : do i = 1, SIZE(omega_in_a)

       n_iter = 0
       n_iter_def = 0

       ! Set up the initial points

       omega_a = c_ext_t(omega_in_a(i))
       omega_b = c_ext_t(omega_in_b(i))
       omega_r = real(0.5_WP*(omega_a + omega_b))

       ! Set up the discriminant function

       st = c_state_t(omega_r=omega_r)
       df = c_discrim_func_t(bp, st, omega_min, omega_max)

       call df%eval(omega_a, discrim_a, status)
       if (status /= STATUS_OK) then
          call report_status_(status, 'initial trial (a)')
          cycle in_loop
       endif
          
       call df%eval(omega_b, discrim_b, status)
       if (status /= STATUS_OK) then
          call report_status_(status, 'initial trial (b)')
          cycle in_loop
       endif
          
       ! If necessary, do a preliminary root find using the deflated
       ! discriminant

       if (nm_p%deflate_roots) then

          st_def = c_state_t(omega_r=omega_r)
          df_def = c_discrim_func_t(bp, st, omega_min, omega_max, omega_def)

          call narrow(df_def, nm_p, omega_a, omega_b, r_ext_t(0._WP), status, n_iter=n_iter_def, n_iter_max=nm_p%n_iter_max)
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate narrow')
             cycle in_loop
          endif

          ! If necessary, reset omega_a and omega_b so they are not
          ! coincident; and then save the revised discriminant values

          if(omega_b == omega_a) then
             omega_b = omega_a*(1._WP + EPSILON(0._WP)*(omega_a/abs(omega_a)))
          endif

          call expand(df, omega_a, omega_b, r_ext_t(0._WP), status, f_cx_a=discrim_a_rev, f_cx_b=discrim_b_rev) 
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate re-expand')
             cycle in_loop
          endif

       else

          discrim_a_rev = discrim_a
          discrim_b_rev = discrim_b

          n_iter_def = 0

       endif

       ! Find the discriminant root

       call solve(df, nm_p, omega_a, omega_b, r_ext_t(0._WP), omega_root, status, &
                  n_iter=n_iter, n_iter_max=nm_p%n_iter_max-n_iter_def, f_cx_a=discrim_a_rev, f_cx_b=discrim_b_rev)
       if (status /= STATUS_OK) then
          call report_status_(status, 'solve')
          cycle in_loop
       endif

       ! Construct the mode_t

       st = c_state_t(omega=cmplx(omega_root), omega_r=omega_r)

       select type (bp)
       type is (nad_bvp_t)
          wv = wave_t(bp, st, j_in(i))
       class default
          $ABORT(Invalid bp class)
       end select

       md = mode_t(wv)

       ! Process it

       chi = abs(md%discrim)/max(abs(discrim_a), abs(discrim_b))
       
       call process_mode(md, n_iter_def+n_iter, chi)

       ! Store the frequency in the deflation array

       omega_def = [omega_def,cmplx(omega_root)]

    end do in_loop

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

         write(OUTPUT_UNIT, 110) 'n_iter_def :', n_iter_def
         write(OUTPUT_UNIT, 110) 'n_iter     :', n_iter
110      format(3X,A,1X,I0)

         write(OUTPUT_UNIT, 120) 'omega_in_a :', cmplx(omega_in_a(i))
         write(OUTPUT_UNIT, 120) 'omega_in_b :', cmplx(omega_in_b(i))
         write(OUTPUT_UNIT, 120) 'omega_a    :', cmplx(omega_a)
         write(OUTPUT_UNIT, 120) 'omega_b    :', cmplx(omega_b)
120      format(3X,A,1X,2E24.16)

      end if

      ! Finish

      return

    end subroutine report_status_

  end subroutine prox_search_2_

end module gyre_prox_search
