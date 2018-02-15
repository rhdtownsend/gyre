! Module   : gyre_r_search
! Purpose  : mode searching (real)
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_r_search

  ! Uses

  use core_kinds
  use core_memory
  use core_order
  use core_parallel

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_context
  use gyre_discrim_func
  use gyre_ext
  use gyre_freq
  use gyre_grid
  use gyre_mode
  use gyre_model
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_bvp
  use gyre_root
  use gyre_rot
  use gyre_rot_factory
  use gyre_scan_par
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

  public :: build_scan
  public :: check_scan
  public :: scan_search

contains

  subroutine build_scan (cx, md_p, os_p, sc_p, omega)

    type(context_t), intent(in)        :: cx
    type(mode_par_t), intent(in)       :: md_p
    type(osc_par_t), intent(in)        :: os_p
    type(scan_par_t), intent(in)       :: sc_p(:)
    real(WP), allocatable, intent(out) :: omega(:)

    integer               :: n_omega
    integer               :: i
    real(WP)              :: omega_min
    real(WP)              :: omega_max
    real(WP)              :: freq_g_min
    real(WP)              :: freq_g_max
    real(WP), allocatable :: freq_g(:)
    integer               :: j

    $ASSERT(SIZE(sc_p) >=1,Empty scan_par_t)

    ! Build the frequency scan

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Building frequency scan'
100    format(A)
    endif

    ! Loop through scan_par_t

    n_omega = 0

    allocate(omega(0))

    sc_p_loop : do i = 1,SIZE(sc_p)

       associate (freq_min => sc_p(i)%freq_min, &
                  freq_max => sc_p(i)%freq_max, &
                  n_freq => sc_p(i)%n_freq, &
                  freq_min_units => sc_p(i)%freq_min_units, &
                  freq_max_units => sc_p(i)%freq_max_units, &
                  freq_frame => sc_p(i)%freq_frame, &
                  grid_frame => sc_p(i)%grid_frame, &
                  grid_type => sc_p(i)%grid_type)
         
         ! Calculate the dimensionless frequency range in the inertial frame
         
         omega_min = omega_from_freq(freq_min, cx%ml, cx%pt_i, cx%pt_o, freq_min_units, freq_frame, md_p, os_p)
         omega_max = omega_from_freq(freq_max, cx%ml, cx%pt_i, cx%pt_o, freq_max_units, freq_frame, md_p, os_p)

         ! Check that the range is valid

         if (omega_max > omega_min) then

            ! Calculate the frequency range in the grid frame

            freq_g_min = freq_from_omega(omega_min, cx%ml, cx%pt_i, cx%pt_o, 'NONE', grid_frame, md_p, os_p)
            freq_g_max = freq_from_omega(omega_max, cx%ml, cx%pt_i, cx%pt_o, 'NONE', grid_frame, md_p, os_p)

            ! Set up the frequencies

            allocate(freq_g(n_freq))

            do j = 1, n_freq
             
               select case(grid_type)
               case('LINEAR')
                  freq_g(j) = ((n_freq-j)*freq_g_min + (j-1)*freq_g_max)/(n_freq-1)
               case('INVERSE')
                  freq_g(j) = (n_freq-1)/((n_freq-j)/freq_g_min + (j-1)/freq_g_max)
               case default
                  $ABORT(Invalid grid_type)
               end select
               
            end do

            ! Store them

            call reallocate(omega, [n_omega+n_freq])

            do j = 1, n_freq
               omega(n_omega+j) = omega_from_freq(freq_g(j), cx%ml, cx%pt_i, cx%pt_o, 'NONE', grid_frame, md_p, os_p)
            end do
          
            n_omega = n_omega + n_freq

            deallocate(freq_g)

            if (check_log_level('INFO')) then
               write(OUTPUT_UNIT, 110) 'added scan interval : ', omega_min, ' -> ', omega_max, ' (', n_freq, ' points, ', TRIM(grid_type), ')'
110            format(3X,A,E11.4,A,E11.4,A,I0,A,A,A)
            endif

         else

            if (check_log_level('INFO')) then
               write(OUTPUT_UNIT, 120) 'ignoring scan interval :', omega_min, ' -> ', omega_max
120            format(3X,A,E24.16,A,E24.16)
            endif

         endif

       end associate

    end do sc_p_loop

    ! Sort the frequencies

    omega = omega(sort_indices(omega))

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine build_scan

  !****

  subroutine check_scan (ml, gr, omega, md_p, os_p)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    real(WP), intent(in)                :: omega(:)
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p

    real(WP) :: Omega_rot
    real(WP) :: omega_c(gr%n_k)
    real(WP) :: omega_c_prev(gr%n_k)
    integer  :: j
    integer  :: k

    ! Check the frequency scan to ensure no zero crossings in omega_c
    ! arise

    if (SIZE(omega) >= 1) then

       !$OMP PARALLEL DO PRIVATE (Omega_rot)
       do k = 1, gr%n_k
          Omega_rot = ml%coeff(I_OMEGA_ROT, gr%pt(k))
          omega_c(k) = omega_corot(omega(1), Omega_rot, md_p%m)
       end do

       $ASSERT(ALL(omega_c > 0._WP) .OR. ALL(omega_c < 0._WP),Critical layer encountered)

       do j = 2, SIZE(omega)

          omega_c_prev = omega_c

          !$OMP PARALLEL DO PRIVATE (Omega_rot)
          do k = 1, gr%n_k
             Omega_rot = ml%coeff(I_OMEGA_ROT, gr%pt(k))
             omega_c(k) = omega_corot(omega(j), Omega_rot, md_p%m)
          end do

          $ASSERT(ALL(SIGN(1._WP, omega_c) == SIGN(1._WP, omega_c_prev)),Transition between prograde and retrograde)

       end do

    endif

    ! Finish

    return

  end subroutine check_scan

  !****

  subroutine scan_search (bp, omega, omega_min, omega_max, process_mode, nm_p)

    class(r_bvp_t), target, intent(inout) :: bp
    real(WP), intent(in)                  :: omega(:)
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
    type(num_par_t), intent(in)           :: nm_p

    real(WP), allocatable      :: omega_a(:)
    real(WP), allocatable      :: omega_b(:)
    type(r_ext_t), allocatable :: discrim_a(:)
    type(r_ext_t), allocatable :: discrim_b(:)
    type(r_state_t)            :: st
    type(r_discrim_func_t)     :: df
    integer                    :: c_beg
    integer                    :: c_end
    integer                    :: c_rate
    integer                    :: n_iter
    type(r_ext_t)              :: omega_root
    integer                    :: status
    integer                    :: i
    type(wave_t)               :: wv
    type(mode_t)               :: md
    type(r_ext_t)              :: chi

    ! Find discriminant root brackets

    call find_root_brackets(bp, omega, omega_min, omega_max, omega_a, omega_b, discrim_a, discrim_b)

    ! Set up the discriminant function

    st = r_state_t(omega=0._WP)
    df = r_discrim_func_t(bp, st, omega_min, omega_max)

    ! Process each bracket to find modes

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) 'Root Solving'
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'm', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(1X,A3,1X,A4,1X,A7,1X,A6,1X,A6,1X,A15,1X,A15,1X,A10,1X,A6)
       
    endif

    call SYSTEM_CLOCK(c_beg, c_rate)

    mode_loop : do i = 1, SIZE(omega_a)

       n_iter = 0

       ! Find the discriminant root

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
          wv = wave_t(bp, st)
       type is (rad_bvp_t)
          wv = wave_t(bp, st)
       class default
          $ABORT(Invalid bp class)
       end select

       md = mode_t(wv, j_m)

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
      
  end subroutine scan_search

  !****

  subroutine find_root_brackets (bp, omega, omega_min, omega_max, omega_a, omega_b, discrim_a, discrim_b)

    class(r_bvp_t), target, intent(inout)   :: bp
    real(WP), intent(in)                    :: omega(:)
    real(WP), intent(in)                    :: omega_min
    real(WP), intent(in)                    :: omega_max
    real(WP), allocatable, intent(out)      :: omega_a(:)
    real(WP), allocatable, intent(out)      :: omega_b(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_a(:)
    type(r_ext_t), allocatable, intent(out) :: discrim_b(:)

    type(r_state_t)        :: st
    type(r_discrim_func_t) :: df
    integer                :: n_omega
    integer                :: c_beg
    integer                :: c_end
    integer                :: c_rate
    integer                :: i
    type(r_ext_t)          :: discrim(SIZE(omega))
    integer                :: status(SIZE(omega))
    integer                :: n_brack
    integer                :: i_brack(SIZE(omega))

    ! Set up the discriminant function

    st = r_state_t(omega=0._WP)
    df = r_discrim_func_t(bp, st, omega_min, omega_max)

    ! Calculate the discriminant on the omega abscissa

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Root bracketing'
100    format(A)
    endif

    n_omega = SIZE(omega)

    call SYSTEM_CLOCK(c_beg, c_rate)

    discrim_loop : do i = 1, n_omega

       call df%eval(r_ext_t(omega(i)), discrim(i), status(i))

       if (check_log_level('DEBUG')) then
          write(OUTPUT_UNIT, 110) omega(i), fraction(discrim(i)), exponent(discrim(i))
110       format(2X,E24.16,2X,F19.16,2X,I7)
       endif

    end do discrim_loop

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
120    format(2X,A,F10.3,1X,A/)
    endif

    ! Find root brackets

    n_brack = 0

    bracket_loop : do i = 1, n_omega-1

       if (status(i) == STATUS_OK .AND. status(i+1) == STATUS_OK) then
          if (discrim(i)*discrim(i+1) <= r_ext_t(0._WP)) then
             n_brack = n_brack + 1
             i_brack(n_brack) = i
          end if
       end if

    end do bracket_loop

    ! Set up the bracket frequencies

    omega_a = omega(i_brack(:n_brack))
    omega_b = omega(i_brack(:n_brack)+1)

    discrim_a = discrim(i_brack(:n_brack))
    discrim_b = discrim(i_brack(:n_brack)+1)

    ! Finish

    return

  end subroutine find_root_brackets

! !****

!   subroutine min_search (bp, omega, process_mode)

!     class(bvp_t), intent(inout) :: bp
!     real(WP), intent(in)        :: omega(:)
!     interface
!        subroutine process_mode (md)
!          use gyre_mode
!          type(mode_t), intent(in) :: md
!        end subroutine process_mode
!     end interface

!     real(WP), allocatable         :: omega_a(:)
!     real(WP), allocatable         :: omega_b(:)
!     real(WP), allocatable         :: omega_c(:)
!     type(ext_real_t), allocatable :: discrim_a(:)
!     type(ext_real_t), allocatable :: discrim_b(:)
!     type(ext_real_t), allocatable :: discrim_c(:)
!     integer                       :: n_brack
!     integer                       :: c_beg
!     integer                       :: c_end
!     integer                       :: c_rate
!     real(WP)                      :: domega
!     type(mode_t)                  :: md
!     integer                       :: i

!     ! Find discriminant minimum brackets

!     call find_min_brackets(bp, omega, omega_a, omega_b, omega_c, discrim_a, discrim_b, discrim_c)

!     ! Process each bracket to find modes

!     if (check_log_level('INFO')) then

!        write(OUTPUT_UNIT, 100) 'Root Solving'
! 100    format(A)

!        write(OUTPUT_UNIT, 110) 'l', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter', 'n'
! 110    format(4(2X,A8),3(2X,A24),2X,A6,2X,A7)
       
!     endif

!     call SYSTEM_CLOCK(c_beg, c_rate)

!     mode_loop : do i = 1, SIZE(omega_a)

!        ! Find the mode

!        domega = MIN(ABS(omega_b(i)-omega_a(i)), ABS(omega_b(i)-omega_c(i)))

!        md = bp%mode([CMPLX(omega_b(i),  0.5_WP*domega, KIND=WP),  &
!                      CMPLX(omega_b(i), -0.5_WP*domega, KIND=WP)])
 
!        if (md%n_pg < md%mp%n_pg_min .OR. md%n_pg > md%mp%n_pg_max) cycle mode_loop

!        ! Process it

!        if (check_log_level('INFO')) then
!           write(OUTPUT_UNIT, 120) md%mp%l, md%n_pg, md%n_p, md%n_g, md%omega, real(md%chi), md%n_iter, md%n
! 120       format(4(2X,I8),3(2X,E24.16),2X,I6,2X,I7)
!        endif

!        call process_mode(md)

!     end do mode_loop

!     call SYSTEM_CLOCK(c_end)

!     if (check_log_level('INFO')) then
!        write(OUTPUT_UNIT, 130) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
! 130    format(2X,A,1X,F10.3,1X,A)
!     endif

!     ! Finish

!     return

!   end subroutine min_search

!****

!   subroutine find_min_brackets (df, omega, omega_a, omega_b, omega_c, discrim_a, discrim_b, discrim_c)

!     type(discrim_func_t), intent(inout)        :: df
!     real(WP), intent(in)                       :: omega(:)
!     real(WP), allocatable, intent(out)         :: omega_a(:)
!     real(WP), allocatable, intent(out)         :: omega_b(:)
!     real(WP), allocatable, intent(out)         :: omega_c(:)
!     type(ext_real_t), allocatable, intent(out) :: discrim_a(:)
!     type(ext_real_t), allocatable, intent(out) :: discrim_b(:)
!     type(ext_real_t), allocatable, intent(out) :: discrim_c(:)

!     integer          :: n_omega
!     integer          :: c_beg
!     integer          :: c_end
!     integer          :: c_rate
!     integer          :: i
!     type(ext_real_t) :: discrim(SIZE(omega))
!     integer          :: n_brack
!     integer          :: i_brack(SIZE(omega))

!     ! Calculate the discriminant on the omega abscissa

!     if (check_log_level('INFO')) then
!        write(OUTPUT_UNIT, 100) 'Root bracketing'
! 100    format(A)
!     endif

!     n_omega = SIZE(omega)

!     call SYSTEM_CLOCK(c_beg, c_rate)

!     discrim_loop : do i = 1, n_omega

!        discrim(i) = abs(df%eval(ext_complex_t(CMPLX(omega(i), KIND=WP))))

!        if (check_log_level('DEBUG')) then
!           write(OUTPUT_UNIT, 110) omega(i), fraction(discrim(i)), exponent(discrim(i))
! 110       format(2X,E24.16,2X,F19.16,2X,I7)
!        endif

!     end do discrim_loop

!     call SYSTEM_CLOCK(c_end)

!     if (check_log_level('INFO')) then
!        write(OUTPUT_UNIT, 120) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
! 120    format(2X,A,F10.3,1X,A/)
!     endif

!     ! Find minimum brackets

!     n_brack = 0

!     bracket_loop : do i = 2, n_omega-1

!        if (discrim(i) < discrim(i-1) .AND. &
!           discrim(i) < discrim(i+1)) then
!           n_brack = n_brack + 1
!           i_brack(n_brack) = i
!        end if

!     end do bracket_loop

!     ! Set up the bracket frequencies

!     omega_a = omega(i_brack(:n_brack)-1)
!     omega_b = omega(i_brack(:n_brack))
!     omega_c = omega(i_brack(:n_brack)+1)

!     discrim_a = discrim(i_brack(:n_brack)-1)
!     discrim_b = discrim(i_brack(:n_brack))
!     discrim_c = discrim(i_brack(:n_brack)+1)

!     ! Finish

!     return

!   end subroutine find_min_brackets

end module gyre_r_search
