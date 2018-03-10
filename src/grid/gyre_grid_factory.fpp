! Module   : gyre_grid_factory
! Purpose  : factory procedures for grid_t type
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

module gyre_grid_factory

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_context
  use gyre_freq
  use gyre_grid
  use gyre_grid_par
  use gyre_grid_spec
  use gyre_grid_util
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_state
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface grid_t
     module procedure grid_t_context_
     module procedure grid_t_grid_spec_
  end interface grid_t

  ! Access specifiers

  private

  public :: grid_t

  ! Procedures

contains

  function grid_t_context_ (cx, omega, gr_p) result (gr)

    type(context_t), pointer, intent(in) :: cx
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(grid_t)                         :: gr

    type(grid_spec_t) :: gs

    ! Construct the grid_t using the supplied context

    gs = grid_spec_t(cx, omega)

    gr = grid_t([gs], gr_p)

    ! Finish

    return

  end function grid_t_context_

  !****

  function grid_t_grid_spec_ (gs, gr_p) result (gr)

    type(grid_spec_t), intent(in) :: gs(:)
    type(grid_par_t), intent(in)  :: gr_p
    type(grid_t)                  :: gr

    class(model_t), pointer :: ml => null()
    integer                 :: i
    integer                 :: s

    $ASSERT_DEBUG(SIZE(gs) >= 1,Too few grid_specs)

    ! Construct the grid_t using the supplied array of
    ! grid_specs. Each grid spec supplies a context and the set of
    ! frequencies to adopt for that context

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Building x grid'
100    format(A)
    endif

    ! Check that the grid_specs are all associated with the same model

    $ASSERT_DEBUG(ASSOCIATED(gs(1)%cx),Null pointer)

    ml => gs(1)%cx%ml

    check_loop : do i = 1, SIZE(gs)
       $ASSERT_DEBUG(ASSOCIATED(gs(i)%cx),Null pointer)
       $ASSERT_DEBUG(ASSOCIATED(gs(i)%cx%ml, ml),Contexts are associated with different models)
    end do check_loop

    ! Create the scaffold grid

    gr = grid_t(ml%grid(), gr_p%x_i, gr_p%x_o)

    ! Add points to the inner region

    call add_inner_(gs, gr_p, gr)

    ! Add points globally

    call add_global_(gs, gr_p, gr)

    ! Report 

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 120) 'Final grid has', gr%s_o()-gr%s_i()+1, 'segment(s) and', gr%n_k, 'point(s):'
120    format(3X,A,1X,I0,1X,A,1X,I0,1X,A)
       
       seg_loop : do s = gr%s_i(), gr%s_o()

          associate( &
               k_i => gr%k_s_i(s), &
               k_o => gr%k_s_o(s))

            write(OUTPUT_UNIT, 130) 'Segment', s, ': x range', gr%pt(k_i)%x, '->', gr%pt(k_o)%x, &
                 '(', k_i, '->', k_o, ')'
130         format(6X,A,1X,I0,1X,A,1X,F6.4,1X,A,1X,F6.4,1X,A,I0,1X,A,1X,I0,A)

          end associate

       end do seg_loop

       write(OUTPUT_UNIT, *)
        
    end if

    ! Finish

  end function grid_t_grid_spec_

  !****

  subroutine add_inner_ (gs, gr_p, gr)

    type(grid_spec_t), intent(in) :: gs(:)
    type(grid_par_t), intent(in)  :: gr_p
    type(grid_t), intent(inout)   :: gr

    integer  :: k_turn_min
    integer  :: k_turn_max
    real(WP) :: x_turn_min
    real(WP) :: x_turn_max
    integer  :: i
    integer  :: j
    integer  :: k_turn
    real(WP) :: x_turn
    real(WP) :: dx
    integer  :: k_max
    integer  :: k
    integer  :: dn(gr%n_k-1)

    ! Add points in the inner region

    ! First, determine the range (in x and cell number) of inner turning
    ! points, across all contexts and frequencies

    k_turn_min = gr%n_k-1
    k_turn_max = 1

    x_turn_min = gr%pt(gr%n_k)%x
    x_turn_max = gr%pt(1)%x

    spec_loop : do i = 1, SIZE(gs)

       associate (cx => gs(i)%cx, &
                  omega => gs(i)%omega)

         !$OMP PARALLEL DO PRIVATE (k_turn, x_turn) REDUCTION (MIN:k_turn_min,x_turn_min) REDUCTION (MAX:k_turn_max,x_turn_max)
         omega_loop : do j = 1, SIZE(omega)

            call find_turn(cx, gr, r_state_t(omega(j)), k_turn, x_turn)

            k_turn_min = MIN(k_turn_min, k_turn)
            k_turn_max = MAX(k_turn_max, k_turn)

            x_turn_min = MIN(x_turn_min, x_turn)
            x_turn_max = MAX(x_turn_max, x_turn)

         end do omega_loop

       end associate

    end do spec_loop

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Found inner turning points, x range', x_turn_min, '->', x_turn_max
110    format(3X,A,1X,F6.4,1X,A,1X,F6.4)
    endif

    ! Now add points (extending from the boundary to the inner turning
    ! point), to ensure that no cell is larger than dx =
    ! x_turn_min/n_inner

    if (gr_p%n_inner > 0) then

       dx = x_turn_min/gr_p%n_inner

       k_max = MIN(k_turn_max, gr%n_k-1)

       !$OMP PARALLEL DO
       cell_loop : do k = 1, k_max

          associate( &
               x_a => gr%pt(k)%x, &
               x_b => gr%pt(k+1)%x)

            dn(k) = FLOOR((x_b - x_a)/dx)

          end associate

       end do cell_loop

       dn(k_max+1:) = 0

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Adding', SUM(dn), 'inner point(s)'
100       format(3X,A,1X,I0,1X,A)
       endif

       ! Add the points

       gr = grid_t(gr, dn)

    endif

    ! Finish

    return

  end subroutine add_inner_

  !****

  subroutine add_global_ (gs, gr_p, gr)

    type(grid_spec_t), intent(in) :: gs(:)
    type(grid_par_t), intent(in)  :: gr_p
    type(grid_t), intent(inout)   :: gr

    integer              :: i_iter
    integer, allocatable :: dn(:)
    integer              :: i
    integer              :: k
    type(point_t)        :: pt
    real(WP)             :: dx

    ! Add points globally 
    
    ! Iterate until no more points need be added
    
    iter_loop : do i_iter = 1, gr_p%n_iter_max

       if (ALLOCATED(dn)) deallocate(dn)
       allocate(dn(gr%n_k-1))

       ! Loop through grid cells

       !$OMP PARALLEL DO PRIVATE (pt, dx, i)
       cell_loop : do k = 1, gr%n_k-1

          associate ( &
               pt_a => gr%pt(k), &
               pt_b => gr%pt(k+1))
                     
            if (pt_a%s == pt_b%s) then

               pt%s = pt_a%s
               pt%x = 0.5_WP*(pt_a%x + pt_b%x)

               ! Loop over grid_specs

               dx = HUGE(0._WP)

               spec_loop : do i = 1, SIZE(gs)

                  associate (cx => gs(i)%cx, &
                             omega => gs(i)%omega)

                    ! Calculate the minimal desired spacing

                    dx = MIN(dx_dispersion_(cx, pt, omega, gr_p, pt_a%x==0), &
                             dx_thermal_(cx, pt, omega, gr_p), &
                             dx_struct_(cx, pt, gr_p), &
                             dx)

                  end associate

               end do spec_loop

               ! Place a floor on dx

               dx = MAX(dx, gr_p%dx_min)

               ! Set dn

               dn(k) = CEILING((pt_b%x - pt_a%x)/dx) - 1

            else

               dn(k) = 0

            endif

          end associate

       end do cell_loop

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Adding', SUM(dn), 'global point(s) in iteration', i_iter
100       format(3X,A,1X,I0,1X,A,1X,I0)
       endif

       ! Check for completion

       if (ALL(dn == 0)) exit iter_loop

       ! Add points

       gr = grid_t(gr, dn)

    end do iter_loop

    ! Finish

    return

  end subroutine add_global_

  !****

  function dx_dispersion_ (cx, pt, omega, gr_p, origin) result (dx)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega(:)
    type(grid_par_t), intent(in) :: gr_p
    logical, intent(in)          :: origin
    real(WP)                     :: dx

    real(WP)        :: V_g
    real(WP)        :: As
    real(WP)        :: U
    real(WP)        :: c_1
    real(WP)        :: Omega_rot
    real(WP)        :: Omega_rot_i
    real(WP)        :: k_r_real
    real(WP)        :: k_r_imag
    integer         :: j
    type(r_state_t) :: st
    real(WP)        :: omega_c
    real(WP)        :: lambda
    real(WP)        :: l_i
    real(WP)        :: g_0
    real(WP)        :: g_2
    real(WP)        :: g_4
    real(WP)        :: gamma
    real(WP)        :: dx_real
    real(WP)        :: dx_imag

    ! Evaluate the target grid spacing dx at point pt from a local
    ! wave dispersion analysis. If k_r is the local radial wavenumber,
    ! then dx = 2pi MIN( 1./(alpha_osc*REAL(k_r)),  1./(alpha_exp*IMAG(k_r)) ].
    ! This corresponds (approximately) to dx being 1/alpha_[osc|exp] times
    ! the [oscillatory|exponential] wavelength
    !
    ! Since k_r depends on the oscillation frequency, it is evaluated
    ! over the supplied range of frequencies omega, and the maximum
    ! real and imaginary parts taken

    if (gr_p%alpha_osc > 0._WP .OR. gr_p%alpha_exp > 0._WP) then

       ! Evaluate coefficients
       
       associate (ml => cx%ml)

         V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)

         Omega_rot = ml%coeff(I_OMEGA_ROT, pt)
         Omega_rot_i = ml%coeff(I_OMEGA_ROT, cx%pt_i)

       end associate

       ! Loop over omega, finding the maximum k_r_real and k_r_imag

       k_r_real = 0._WP
       k_r_imag = 0._WP

!       !$OMP PARALLEL DO PRIVATE (st, omega_c, lambda, l_i, g_0, g_2, g_4, gamma) REDUCTION (MAX:k_r_real,k_r_imag)
       omega_loop : do j = 1, SIZE(omega)

          st = r_state_t(omega(j))

          omega_c = cx%omega_c(Omega_rot, st)

          lambda = cx%lambda(Omega_rot, st)
          l_i = cx%l_e(Omega_rot_i, st)

          ! Calculate the propagation discriminant gamma

          g_4 = -4._WP*V_g*c_1
          g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
          g_0 = -4._WP*lambda*As/c_1

          gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

          ! Update the maximal k_r

          if (gamma < 0._WP) then

             ! Propagation zone

             k_r_real = MAX(k_r_real, ABS(0.5_WP*SQRT(-gamma))/pt%x)
             k_r_imag = MAX(k_r_imag, ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i))/pt%x)

          else

             ! Evanescent zone; if we're adjacent to the origin, drop the divering root

             if (origin) then
                k_r_imag = MAX(k_r_imag, ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i + SQRT(gamma)))/pt%x)
             else
                k_r_imag = MAX(k_r_imag, ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i - SQRT(gamma)))/pt%x, & 
                                         ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i + SQRT(gamma)))/pt%x)
             endif

          end if

       end do omega_loop

       ! Now calculate dx

       if (k_r_real > 0._WP) then
          dx_real = TWOPI/(gr_p%alpha_osc*k_r_real)
       else
          dx_real = HUGE(0._WP)
       endif

       if (k_r_imag > 0._WP) then
          dx_imag = TWOPI/(gr_p%alpha_osc*k_r_imag)
       else
          dx_imag = HUGE(0._WP)
       endif

       dx = MIN(dx_real, dx_imag)

    else

       dx = HUGE(0._WP)

    end if

    ! Finish

    return

  end function dx_dispersion_

  !****

  function dx_thermal_ (cx, pt, omega, gr_p) result (dx)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega(:)
    type(grid_par_t), intent(in) :: gr_p
    real(WP)                     :: dx

    real(WP)        :: V
    real(WP)        :: nabla
    real(WP)        :: c_rad
    real(WP)        :: c_thk
    real(WP)        :: Omega_rot
    real(WP)        :: tau
    integer         :: j
    type(r_state_t) :: st
    real(WP)        :: omega_c

    ! Evaluate the target grid spacing dx at point pt from a local
    ! thermal dispersion analysis. If tau is the (absolute) value of
    ! the eigenvalue associated with thermal parts of the
    ! non-adiabatic equations, wavenumber, then dx = 2pi/(alpha_thm*tau)
    !
    ! Since tau depends on the oscillation frequency, it is evaluated
    ! over the supplied range of frequencies omega, and the maximum
    ! absolute value taken

    if (gr_p%alpha_thm > 0._WP) then

       ! Evaluate coefficients

       associate (ml => cx%ml)

         V = ml%coeff(I_V_2, pt)*pt%x**2
         nabla = ml%coeff(I_NABLA, pt)

         c_rad = ml%coeff(I_C_RAD, pt)
         c_thk = ml%coeff(I_C_THK, pt)

         Omega_rot = ml%coeff(I_OMEGA_ROT, pt)

       end associate

       ! Loop over omega, finding the maximum tau

       tau = 0._WP

!       !$OMP PARALLEL DO PRIVATE (st, omega_c) REDUCTION (MAX:tau)
       omega_loop : do j = 1, SIZE(omega)

          st = r_state_t(omega(j))

          omega_c = cx%omega_c(Omega_rot, st)
         
          ! Update the maximal tau
          
          tau = MAX(tau, SQRT(ABS(V*nabla*omega_c*c_thk/c_rad))/pt%x)

       end do omega_loop

       ! Now calculate dx

       dx = TWOPI/(gr_p%alpha_thm*tau)

    else

       dx = HUGE(0._WP)

    end if

    ! Finish

    return

  end function dx_thermal_

  !****

  function dx_struct_ (cx, pt, gr_p) result (dx)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt
    type(grid_par_t), intent(in) :: gr_p
    real(WP)                     :: dx

    real(WP) :: dV_2
    real(WP) :: dAs
    real(WP) :: dU
    real(WP) :: dc_1
    real(WP) :: dGamma_1

    ! Evaluate the target grid spacing dx at sample points pt to
    ! ensure adequate model structure resolution. For a single
    ! structure coefficient C, dx = x/(alpha_str*dlnC/dlnx); we take
    ! the minimum dx over a number of C

    dx = HUGE(0._WP)

    if (gr_p%alpha_str > 0._WP) then

       ! Evaluate coefficients

       associate (ml => cx%ml)

         dV_2 = ml%dcoeff(I_V_2, pt)
         dAs = ml%dcoeff(I_AS, pt)
         dU = ml%dcoeff(I_U, pt)
         dc_1 = ml%dcoeff(I_C_1, pt)
         dGamma_1 = ml%dcoeff(I_GAMMA_1, pt)

       end associate

       if (dV_2 /= 0._WP) dx = MIN(dx, ABS(pt%x/(gr_p%alpha_str*dV_2)))
       if (dAs /= 0._WP) dx = MIN(dx, ABS(pt%x/(gr_p%alpha_str*dAs)))
       if (dU /= 0._WP) dx = MIN(dx, ABS(pt%x/(gr_p%alpha_str*dU)))
       if (dc_1 /= 0._WP) dx = MIN(dx, ABS(pt%x/(gr_p%alpha_str*dc_1)))
       if (dGamma_1 /= 0._WP) dx = MIN(dx, ABS(pt%x/(gr_p%alpha_str*dGamma_1)))

    else

       dx = HUGE(0._WP)

    end if

    ! Finish

    return

  end function dx_struct_

end module gyre_grid_factory
