! Module  : prox_search_m
! Purpose : mode searching (complex, proximity)
!
! Copyright 2013-2022 Rich Townsend & The GYRE Team
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

module prox_search_m

  ! Uses

  use kinds_m
  use constants_m

  use bvp_m
  use discrim_m
  use ext_m
  use math_m
  use mode_m
  use nad_bvp_m
  use num_par_m
  use root_m
  use state_m
  use status_m
  use tnad_bvp_m
  use util_m
  use wave_m

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

  subroutine prox_search_1_ (bp, omega_in, id_in, omega_min, omega_max, nm_p, process_mode)

    class(c_bvp_t), target, intent(inout) :: bp
    complex(WP), intent(in)               :: omega_in(:)
    integer, intent(in)                   :: id_in(:)
    real(WP), intent(in)                  :: omega_min
    real(WP), intent(in)                  :: omega_max
    type(num_par_t), intent(in)           :: nm_p
    interface
       subroutine process_mode (md, n_iter, chi)
         use kinds_m
         use ext_m
         use mode_m
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

    $CHECK_BOUNDS(SIZE(id_in),SIZE(omega_in))

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
    
    call prox_search_2_(bp, omega_in_a, omega_in_b, id_in, omega_min, omega_max, nm_p, process_mode)

    ! Finish

    return

  end subroutine prox_search_1_

  !****

  subroutine prox_search_2_ (bp, omega_in_a, omega_in_b, id_in, omega_min, omega_max, nm_p, process_mode)

    class(c_bvp_t), intent(inout) :: bp
    complex(WP), intent(in)       :: omega_in_a(:)
    complex(WP), intent(in)       :: omega_in_b(:)
    integer, intent(in)           :: id_in(:)
    real(WP), intent(in)          :: omega_min
    real(WP), intent(in)          :: omega_max
    type(num_par_t), intent(in)   :: nm_p
    interface
       subroutine process_mode (md, n_iter, chi)
         use kinds_m
         use ext_m
         use mode_m
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
    integer                  :: n_iter
    integer                  :: n_iter_def
    real(WP)                 :: omega_r
    type(c_ext_t)            :: discrim_a
    type(c_ext_t)            :: discrim_b
    type(c_ext_t)            :: z_a
    type(c_ext_t)            :: z_b
    type(c_ext_t)            :: f_z_a
    type(c_ext_t)            :: f_z_b
    logical                  :: apply_deflate
    integer                  :: status
    type(c_ext_t)            :: z_root
    type(wave_t)             :: wv
    type(mode_t)             :: md
    type(r_ext_t)            :: chi

    $CHECK_BOUNDS(SIZE(id_in),SIZE(omega_in_a))
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

       omega_r = real(0.5_WP*(omega_in_a(i) + omega_in_b(i)))

       call eval_discrim(bp, c_state_t(omega_in_a(i), omega_r=omega_r), discrim_a)
       call eval_discrim(bp, c_state_t(omega_in_b(i), omega_r=omega_r), discrim_b)

       z_a = c_ext_t(omega_in_a(i))
       z_b = c_ext_t(omega_in_b(i))
       f_z_a = discrim_a
       f_z_b = discrim_b

       ! If necessary, do a preliminary root find using the deflated
       ! discriminant

       if (nm_p%deflate_roots) then

          apply_deflate = .TRUE.

          call narrow_bracket(eval_discrim_z_, z_a, z_b, r_ext_t(0._WP), nm_p, &
               status, n_iter=n_iter_def, n_iter_max=nm_p%n_iter_max)
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate narrow')
             cycle in_loop
          endif

          ! If necessary, reset z_a and z_b so they are not
          ! coincident; and then save the revised discriminant values

          if (z_b == z_a) then
             z_b = z_a*(1._WP + EPSILON(0._WP)*(z_a/abs(z_a)))
          endif

          call expand_bracket(eval_discrim_z_, z_a, z_b, r_ext_t(0._WP), &
               status, f_z_a=f_z_a, f_z_b=f_z_b)
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate re-expand')
             cycle in_loop
          endif

       else

          n_iter_def = 0

       endif

       ! Find the discriminant root

       apply_deflate = .FALSE.

       call solve_root(eval_discrim_z_, z_a, z_b, r_ext_t(0._WP), nm_p, &
            z_root, status, n_iter=n_iter, n_iter_max=nm_p%n_iter_max-n_iter_def, &
            f_z_a=f_z_a, f_z_b=f_z_b)
       if (status /= STATUS_OK) then
          call report_status_(status, 'solve')
          cycle in_loop
       endif

       ! Construct the mode_t

       select type (bp)
       type is (nad_bvp_t)
          wv = wave_t(bp, c_state_t(cmplx(z_root), omega_r=omega_r), id_in(i))
       type is (tnad_bvp_t)
          wv = wave_t(bp, c_state_t(cmplx(z_root), omega_r=omega_r), id_in(i))
       class default
          $ABORT(Invalid bp class)
       end select

       md = mode_t(wv)

       ! Process it

       chi = abs(md%discrim)/max(abs(discrim_a), abs(discrim_b))
       
       call process_mode(md, n_iter_def+n_iter, chi)

       ! Store the frequency in the deflation array

       omega_def = [omega_def,cmplx(z_root)]

    end do in_loop

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(2X,A,1X,F10.3,1X,A/)
    endif

    ! Finish

    return

  contains

    subroutine eval_discrim_z_ (z, discrim, status)

      type(c_ext_t), intent(in)  :: z
      type(c_ext_t), intent(out) :: discrim
      integer, intent(out)       :: status

      complex(WP) :: omega

      ! Evaluate the frequency from z

      omega = cmplx(z)

      ! Evaluate the discriminant

      if (REAL(omega) >= omega_min .AND. REAL(omega) <= omega_max) then

         call eval_discrim(bp, c_state_t(omega, omega_r=omega_r), discrim)

         if (apply_deflate) then
            discrim = discrim*PRODUCT(omega/(omega - omega_def))
         end if
            
         status = STATUS_OK

      else

         status = STATUS_OMEGA_DOMAIN

      endif
      
      ! Finish

      return

    end subroutine eval_discrim_z_

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

         write(OUTPUT_UNIT, 110) 'n_iter_def :', n_iter_def
         write(OUTPUT_UNIT, 110) 'n_iter     :', n_iter
110      format(3X,A,1X,I0)

         write(OUTPUT_UNIT, 120) 'omega_in_a :', cmplx(omega_in_a(i))
         write(OUTPUT_UNIT, 120) 'omega_in_b :', cmplx(omega_in_b(i))
120      format(3X,A,1X,2E24.16)

      end if

      ! Finish

      return

    end subroutine report_status_

  end subroutine prox_search_2_

end module prox_search_m
