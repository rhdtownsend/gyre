! Module   : gyre_scan
! Purpose  : frequency scanning routines
!
! Copyright 2013-2018 Rich Townsend
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

module gyre_scan

  ! Uses

  use core_kinds
  use core_order

  use gyre_context
  use gyre_freq
  use gyre_grid
  use gyre_mode
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_scan_par
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_scan
  public :: check_scan

contains

  subroutine build_scan (cx, md_p, os_p, sc_p, omega)

    type(context_t), intent(in)        :: cx
    type(mode_par_t), intent(in)       :: md_p
    type(osc_par_t), intent(in)        :: os_p
    type(scan_par_t), intent(in)       :: sc_p(:)
    real(WP), allocatable, intent(out) :: omega(:)

    class(model_t), pointer :: ml
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o
    integer                 :: i
    real(WP)                :: omega_i_min
    real(WP)                :: omega_i_max
    real(WP)                :: omega_g_min
    real(WP)                :: omega_g_max
    real(WP)                :: omega_g
    real(WP), allocatable   :: omega_i(:)
    integer                 :: j

    $ASSERT(SIZE(sc_p) >=1,Empty scan_par_t)

    ! Build the frequency scan

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Building frequency scan'
100    format(A)
    endif

    ! Loop through scan_par_t

    ml => cx%model()
    pt_i = cx%point_i()
    pt_o = cx%point_o()
         
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

         omega_i_min = omega_from_freq(freq_min, ml, pt_i, pt_o, freq_min_units, freq_frame, md_p, os_p)
         omega_i_max = omega_from_freq(freq_max, ml ,pt_i, pt_o, freq_max_units, freq_frame, md_p, os_p)

         ! Check that the range is valid

         if (omega_i_max > omega_i_min) then

            ! Calculate the dimensionless frequency range in the grid frame

            omega_g_min = freq_from_omega(omega_i_min, ml, pt_i, pt_o, 'NONE', grid_frame, md_p, os_p)
            omega_g_max = freq_from_omega(omega_i_max, ml, pt_i, pt_o, 'NONE', grid_frame, md_p, os_p)

            ! Set up the frequencies

            allocate(omega_i(n_freq))

            do j = 1, n_freq

               ! Grid frame

               select case(grid_type)
               case('LINEAR')
                  omega_g = ((n_freq-j)*omega_g_min + (j-1)*omega_g_max)/(n_freq-1)
               case('INVERSE')
                  omega_g = (n_freq-1)/((n_freq-j)/omega_g_min + (j-1)/omega_g_max)
               case default
                  $ABORT(Invalid grid_type)
               end select

               ! Inertial frame

               omega_i(j) = omega_from_freq(omega_g, ml, pt_i, pt_o, 'NONE', grid_frame, md_p, os_p)

            end do

            ! Store them

            omega = [omega,omega_i]

            if (check_log_level('INFO')) then
               write(OUTPUT_UNIT, 110) 'added scan interval : ', omega_i_min, ' -> ', omega_i_max, ' (', n_freq, ' points, ', TRIM(grid_type), ')'
110            format(3X,A,E11.4,A,E11.4,A,I0,A,A,A)
            endif

            deallocate(omega_i)

         else

            if (check_log_level('INFO')) then
               write(OUTPUT_UNIT, 120) 'ignoring scan interval :', omega_i_min, ' -> ', omega_i_max
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
    real(WP), intent(inout)             :: omega(:)
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p

    real(WP) :: Omega_rot
    real(WP) :: omega_c(gr%n_k)
    real(WP) :: omega_c_prev(gr%n_k)
    integer  :: j
    integer  :: k

    ! Check the frequency scan to remove any points with leaky
    ! boundary conditions

    ! This needs to be done

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

end module gyre_scan
