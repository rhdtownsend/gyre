! Module   : gyre_c_search
! Purpose  : mode searching (complex)
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_c_search

  ! Uses

  use core_kinds
  use gyre_constants
  use core_order
  use core_parallel

  use gyre_bvp
  use gyre_discrim_func
  use gyre_ext
  use gyre_mode
  use gyre_mode_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_root
  use gyre_state
  use gyre_status
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: prox_search

contains

  subroutine prox_search (bp, md_in, omega_min, omega_max, process_mode, md_p, nm_p, os_p)

    class(c_bvp_t), target, intent(inout) :: bp
    type(mode_t), intent(in)              :: md_in(:)
    real(WP), intent(in)                  :: omega_min
    real(WP), intent(in)                  :: omega_max
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
    type(mode_par_t), intent(in)          :: md_p
    type(num_par_t), intent(in)           :: nm_p
    type(osc_par_t), intent(in)           :: os_p
    
    complex(WP), allocatable :: omega_def(:)
    integer                  :: c_beg
    integer                  :: c_end
    integer                  :: c_rate
    integer                  :: n_md_in
    integer                  :: i
    real(WP)                 :: domega
    type(c_ext_t)            :: omega_a
    type(c_ext_t)            :: omega_b
    real(WP)                 :: omega_ad
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

    ! Initialize the frequency deflation array

    allocate(omega_def(0))

    ! Process each initial mode to find a proximate mode

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'm', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(1X,A3,1X,A4,1X,A7,1X,A6,1X,A6,1X,A15,1X,A15,1X,A10,1X,A6)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    n_md_in = SIZE(md_in)

    mode_loop : do i = 1, n_md_in

       n_iter = 0
       n_iter_def = 0

       ! Set up initial guesses

       if (n_md_in > 1) then

          if (i == 1) then
             domega = ABS(md_in(2)%omega - md_in(1)%omega)
          elseif (i == n_md_in) then
             domega = ABS(md_in(n_md_in)%omega - md_in(n_md_in-1)%omega)
          else
             domega = MIN(ABS(md_in(i)%omega - md_in(i-1)%omega), &
                          ABS(md_in(i+1)%omega - md_in(i)%omega))
          endif

          domega = domega*1E-3 ! This seems rather random

       else

          domega = ABS(md_in(i)%omega)*SQRT(EPSILON(0._WP))

       endif

       omega_a = c_ext_t(md_in(i)%omega + CMPLX(0._WP, domega, KIND=WP))
       omega_b = c_ext_t(md_in(i)%omega - CMPLX(0._WP, domega, KIND=WP))

       omega_ad = real(md_in(i)%omega)

       ! call improve_omega(bp, md_p, os_p, md_in(i)%x, omega_a)
       ! call improve_omega(bp, md_p, os_p, md_in(i)%x, omega_b)

       st = c_state_t(omega=0._WP, omega_r=omega_ad)
       df = c_discrim_func_t(bp, st, omega_min, omega_max)

       call df%eval(omega_a, discrim_a, status)
       if (status /= STATUS_OK) then
          call report_status_(status, 'initial guess (a)')
          cycle mode_loop
       endif
          
       call df%eval(omega_b, discrim_b, status)
       if (status /= STATUS_OK) then
          call report_status_(status, 'initial guess (b)')
          cycle mode_loop
       endif
          
       ! If necessary, do a preliminary root find using the deflated
       ! discriminant

       if (nm_p%deflate_roots) then

          st_def = c_state_t(omega=0._WP, omega_r=omega_ad)
          df_def = c_discrim_func_t(bp, st, omega_min, omega_max, omega_def)

          call narrow(df_def, nm_p, omega_a, omega_b, r_ext_t(0._WP), status, n_iter=n_iter_def, n_iter_max=nm_p%n_iter_max)
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate narrow')
             cycle mode_loop
          endif

          ! If necessary, reset omega_a and omega_b so they are not
          ! coincident; and then save the revised discriminant values

          if(omega_b == omega_a) then
             omega_b = omega_a*(1._WP + EPSILON(0._WP)*(omega_a/ABS(omega_a)))
          endif

          call expand(df, omega_a, omega_b, r_ext_t(0._WP), status, f_cx_a=discrim_a_rev, f_cx_b=discrim_b_rev) 
          if (status /= STATUS_OK) then
             call report_status_(status, 'deflate re-expand')
             cycle mode_loop
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
          cycle mode_loop
       endif

       ! Construct the mode_t

       st = c_state_t(omega=cmplx(omega_root), omega_r=omega_ad)

       select type (bp)
       type is (nad_bvp_t)
          wv = wave_t(bp, st)
       class default
          $ABORT(Invalid bp class)
       end select

       md = mode_t(wv, md_in(i)%j)

       ! Process it

       chi = abs(md%discrim)/max(abs(discrim_a), abs(discrim_b))
       
       call process_mode(md, n_iter_def+n_iter, chi)

       ! Store the frequency in the deflation array

       omega_def = [omega_def,cmplx(omega_root)]

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

         write(OUTPUT_UNIT, 110) 'n_iter_def :', n_iter_def
         write(OUTPUT_UNIT, 110) 'n_iter     :', n_iter
110      format(3X,A,1X,I0)

         write(OUTPUT_UNIT, 120) 'omega_a    :', cmplx(omega_a)
         write(OUTPUT_UNIT, 120) 'omega_b    :', cmplx(omega_b)
120      format(3X,A,1X,2E24.16)

      end if

      ! Finish

      return

    end subroutine report_status_

  end subroutine prox_search

!****

  ! subroutine improve_omega (bp, mp, op, x, omega)

  !   class(c_bvp_t), target, intent(inout) :: bp
  !   type(mode_par_t), intent(in)          :: mp
  !   type(osc_par_t), intent(in)           :: op
  !   real(WP), intent(in)                  :: x(:)
  !   type(c_ext_t), intent(inout)          :: omega

  !   integer       :: n
  !   real(WP)      :: x_ref
  !   complex(WP)   :: y(6,SIZE(x))
  !   complex(WP)   :: y_ref(6)
  !   type(c_ext_t) :: discrim
  !   type(mode_t)  :: md

  !   ! Use the integral expression for the eigenfrequency to improve
  !   ! omega

  !   ! Reconstruct on the supplied grid

  !   n = SIZE(x)

  !   x_ref = x(n)

  !   call bp%recon(cmplx(omega), x, x_ref, y, y_ref, discrim)

  !   ! Create the mode

  !   md = mode_t(bp%ml, mp, op, cmplx(omega), discrim, &
  !               x, y, x_ref, y_ref)

  !   ! Improve omega

  !   omega = c_ext_t(md%omega_int())

  !   ! Finish

  !   return

  ! end subroutine improve_omega

end module gyre_c_search
