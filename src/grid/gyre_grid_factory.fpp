! Module   : gyre_grid_factory
! Purpose  : factory procedures for grid_t type
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

module gyre_grid_factory

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_grid
  use gyre_grid_par
  use gyre_grid_util
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface grid_t
     module procedure grid_t_model_
     module procedure grid_t_weights_
  end interface grid_t

  ! Access specifiers

  private

  public :: grid_t

  ! Procedures

contains

  function grid_t_model_ (ml, omega, gr_p, md_p, os_p) result (gr)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(grid_t)                        :: gr

    ! Construct the grid_t using the supplied model grid as the base

    ! Create the scaffold grid

    gr = grid_t(ml%grid(), gr_p%x_i, gr_p%x_o)

    ! Add points at the center

    call add_center_(ml, omega, gr_p, md_p, os_p, gr)

    ! Add points globally

    call add_global_(ml, omega, gr_p, md_p, os_p, gr)

    ! Finish

    return

  end function grid_t_model_

  !****

  function grid_t_weights_ (w, x_i, x_o) result (gr)

    real(WP), intent(in) :: w(:)
    real(WP), intent(in) :: x_i
    real(WP), intent(in) :: x_o
    type(grid_t)         :: gr

    ! Construct the grid_t using the supplied weights array and range
    ! (x_i,x_o)

    gr = grid_t((1._WP-w)*x_i + w*x_o)

    ! Finish

    return

  end function grid_t_weights_

  !****

  subroutine add_center_ (ml, omega, gr_p, md_p, os_p, gr)

    class(model_t), pointer, intent(in)  :: ml
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(grid_t), intent(inout)          :: gr

    integer       :: k_turn
    real(WP)      :: x_turn
    integer       :: j
    integer       :: k_turn_omega
    real(WP)      :: x_turn_omega
    real(WP)      :: dx_max
    integer       :: k
    integer       :: dn(gr%n_k-1)
    type(point_t) :: pt_a
    type(point_t) :: pt_b

    ! Add points at the center of grid gr, to ensure that no cell is
    ! larger than dx_max = x_turn/n_center

    x_turn = HUGE(0._WP)
    k_turn = gr%n_k

    if (gr_p%n_center > 0) then

       ! First, determine the inner turning point (over all
       ! frequencies) closest to the center

       omega_loop : do j = 1, SIZE(omega)

          call find_turn(ml, gr, omega(j), md_p, os_p, k_turn_omega, x_turn_omega)

          if (x_turn_omega < x_turn) then
             k_turn = k_turn_omega
             x_turn = x_turn_omega
          endif

       end do omega_loop

       k_turn = MIN(k_turn, gr%n_k-1)

       ! Add points to the cell containing the turning point, and each
       ! cell inside it, so that none is larger than dx_max

       dx_max = x_turn/gr_p%n_center

       !$OMP PARALLEL DO PRIVATE (pt_a, pt_b)
       cell_loop : do k = 1, k_turn

          pt_a = gr%pt(k)
          pt_b = gr%pt(k+1)
          
          dn(k) = FLOOR((pt_b%x - pt_a%x)/dx_max)

       end do cell_loop

       dn(k_turn+1:) = 0

       gr = grid_t(gr, dn)

    endif

    ! Finish

    return

  end subroutine add_center_

  !****

  subroutine add_global_ (ml, omega, gr_p, md_p, os_p, gr)

    class(model_t), pointer, intent(in)  :: ml
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(grid_t), intent(inout)          :: gr

    integer :: dn(gr%n_k-1)

    ! Add points globally, as determined by the various
    ! grid-resampling parameters in gr_p
    
    dn = MAX(dn_dispersion_(ml, gr, omega, gr_p, md_p, os_p), &
             dn_thermal_(ml, gr, omega, gr_p, md_p, os_p), &
             dn_struct_(ml, gr, gr_p))

    gr = grid_t(gr, dn)

    ! Finish

    return

  end subroutine add_global_

  !****

  function dn_dispersion_ (ml, gr, omega, gr_p, md_p, os_p) result (dn)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer                             :: dn(gr%n_k-1)

    class(r_rot_t), allocatable :: rt
    real(WP)                    :: beta_r_max(gr%n_k)
    real(WP)                    :: beta_i_max(gr%n_k)
    integer                     :: k
    type(point_t)               :: pt
    real(WP)                    :: V_g
    real(WP)                    :: As
    real(WP)                    :: U
    real(WP)                    :: c_1
    integer                     :: j
    real(WP)                    :: omega_c
    real(WP)                    :: lambda
    real(WP)                    :: l_i
    real(WP)                    :: g_4
    real(WP)                    :: g_2
    real(WP)                    :: g_0
    real(WP)                    :: gamma
    type(point_t)               :: pt_a
    type(point_t)               :: pt_b
    real(WP)                    :: dphi_osc
    real(WP)                    :: dphi_exp

    ! Determine how many points dn to add to each cell of the grid,
    ! such there are at least alpha_osc points per oscillatory
    ! wavelength and alpha_exp points per exponential wavelength
    !
    ! Wavelengths are calculated based on a local dispersion analysis
    ! of the adibatic/Cowling wave equation, for inertial frequencies
    ! specified by omega

    if (gr_p%alpha_osc > 0._WP .OR. gr_p%alpha_exp > 0._WP) then

       allocate(rt, SOURCE=r_rot_t(ml, gr, md_p, os_p))

       ! At each point, determine the maximum absolute value of the
       ! real and imaginary parts of the local radial wavenumber beta,
       ! for all possible omega values

       beta_r_max(1) = 0._WP
       beta_i_max(1) = 0._WP

       !$OMP PARALLEL DO PRIVATE (pt, V_g, As, U, c_1, j, omega_c, lambda, l_i, g_4, g_2, g_0, gamma)
       wavenumber_loop : do k = 2, gr%n_k-1

          pt = gr%pt(k)

          V_g = ml%V_2(pt)*pt%x**2/ml%Gamma_1(pt)
          As = ml%As(pt)
          U = ml%U(pt)
          c_1 = ml%c_1(pt)

          beta_r_max(k) = 0._WP
          beta_i_max(k) = 0._WP

          omega_loop : do j = 1, SIZE(omega)

             omega_c = rt%omega_c(pt, omega(j))

             lambda = rt%lambda(pt, omega(j))
             l_i = rt%l_i(omega(j))
            
             ! Calculate the propagation discriminant gamma

             g_4 = -4._WP*V_g*c_1
             g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
             g_0 = -4._WP*lambda*As/c_1

             gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

             ! Update the wavenumber maxima

             if (gamma < 0._WP) then
               
                ! Propagation zone

                beta_r_max(k) = MAX(beta_r_max(k), ABS(0.5_WP*SQRT(-gamma))/pt%x)
                beta_i_max(k) = MAX(beta_i_max(k), ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i))/pt%x)

             else

                ! Evanescent zone

                beta_i_max(k) = MAX(beta_i_max(k), &
                     ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i - SQRT(gamma)))/pt%x, &
                     ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i + SQRT(gamma)))/pt%x)
             
             end if

          end do omega_loop

       end do wavenumber_loop

       beta_r_max(gr%n_k) = 0._WP
       beta_i_max(gr%n_k) = 0._WP

       ! Set up dn

       !$OMP PARALLEL DO PRIVATE (pt_a, pt_b, dphi_osc, dphi_exp)
       cell_loop : do k = 1, gr%n_k-1

          pt_a = gr%pt(k)
          pt_b = gr%pt(k+1)

          if (pt_a%s == pt_b%s) then

             ! Calculate the oscillatory and exponential phase change across
             ! the cell

             dphi_osc = MAX(beta_r_max(k), beta_r_max(k+1))*(pt_b%x - pt_a%x)
             dphi_exp = MAX(beta_i_max(k), beta_i_max(k+1))*(pt_b%x - pt_a%x)

             ! Set up dn

             dn(k) = MAX(FLOOR((gr_p%alpha_osc*dphi_osc)/TWOPI), FLOOR((gr_p%alpha_exp*dphi_exp)/TWOPI))

          else

             dn(k) = 0

          endif

       end do cell_loop

    else

       dn = 0

    endif

    ! Finish

    return

  end function dn_dispersion_

  !****

  function dn_thermal_ (ml, gr, omega, gr_p, md_p, os_p) result (dn)

    class(model_t), pointer, intent(in)  :: ml
    type(grid_t), intent(in)             :: gr
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer                              :: dn(gr%n_k-1)

    type(point_t)               :: pt
    class(r_rot_t), allocatable :: rt
    real(WP)                    :: beta_t_max(gr%n_k)
    integer                     :: k
    real(WP)                    :: V
    real(WP)                    :: nabla
    real(WP)                    :: c_rad
    real(WP)                    :: c_thm
    integer                     :: j
    real(WP)                    :: omega_c
    type(point_t)               :: pt_a
    type(point_t)               :: pt_b
    real(WP)                    :: dphi_thm

    ! Determine how many points dn to add to each cell of the grid,
    ! such there are at least alpha_thm points per thermal length

    if (gr_p%alpha_thm > 0._WP) then

       allocate(rt, SOURCE=r_rot_t(ml, gr, md_p, os_p))

       ! At each point, determine the maximum absolute value of the
       ! local thermal wavenumber beta_t, for all possible omega
       ! values

       beta_t_max(1) = 0._WP

       !$OMP PARALLEL DO PRIVATE (pt, V, nabla, c_rad, c_thm, j, omega_c)
       wavenumber_loop : do k = 2, gr%n_k-1

          pt = gr%pt(k)

          V = ml%V_2(pt)*pt%x**2
          nabla = ml%nabla(pt)

          c_rad = ml%c_rad(pt)
          c_thm = ml%c_thm(pt)

          beta_t_max(k) = 0._WP

          omega_loop : do j = 1, SIZE(omega)

             omega_c = rt%omega_c(pt, omega(j))

             beta_t_max(k) = MAX(beta_t_max(k), SQRT(ABS(V*nabla*omega_c*c_thm/c_rad))/pt%x)

          end do omega_loop

       end do wavenumber_loop

       beta_t_max(gr%n_k) = 0._WP

       ! Set up dn

       !$OMP PARALLEL DO PRIVATE (pt_a, pt_b, dphi_thm)
       cell_loop : do k = 1, gr%n_k-1

          pt_a = gr%pt(k)
          pt_b = gr%pt(k+1)

          if (pt_a%s == pt_b%s) then

             ! Calculate the thermal phase change across the cell

             dphi_thm = MAX(beta_t_max(k), beta_t_max(k+1))*(pt_b%x - pt_a%x)

             ! Set up dn

             dn(k) = FLOOR((gr_p%alpha_thm*dphi_thm)/TWOPI)

          else

             dn(k) = 0

          endif

       end do cell_loop

    else

       dn = 0

    endif

    ! Finish

    return

    return

  end function dn_thermal_

  !****

  function dn_struct_ (ml, gr, gr_p) result (dn)

    class(model_t), pointer, intent(in)  :: ml
    type(grid_t), intent(in)             :: gr
    type(grid_par_t), intent(in)         :: gr_p
    integer                              :: dn(gr%n_k-1)

    integer       :: k
    type(point_t) :: pt_a
    type(point_t) :: pt_b

    ! Determine how many points dn to add to each cell of the grid x,
    ! such there are at least alpha_str points per dex change in the
    ! structure variables (V, As, Gamma_1, c_1, & U). Segment
    ! boundaries are excluded

    if (gr_p%alpha_str > 0) then

       ! Set up dn

       !$OMP PARALLEL DO PRIVATE (pt_a, pt_b)
       cell_loop : do k = 1, gr%n_k-1

          pt_a = gr%pt(k)
          pt_b = gr%pt(k+1)

          if (pt_a%s == pt_b%s) then

             dn(k) = 0

             dn(k) = MAX(dn(k), FLOOR(gr_p%alpha_str*dlog_(ml%V_2(pt_a), ml%V_2(pt_b))))
             dn(k) = MAX(dn(k), FLOOR(gr_p%alpha_str*dlog_(ml%As(pt_a), ml%As(pt_b))))
             dn(k) = MAX(dn(k), FLOOR(gr_p%alpha_str*dlog_(ml%Gamma_1(pt_a), ml%Gamma_1(pt_b))))
             dn(k) = MAX(dn(k), FLOOR(gr_p%alpha_str*dlog_(ml%c_1(pt_a), ml%c_1(pt_b))))
             dn(k) = MAX(dn(k), FLOOR(gr_p%alpha_str*dlog_(ml%U(pt_a), ml%U(pt_b))))

          else

             dn(k) = 0

          endif

       end do cell_loop

    else

       dn = 0

    end if

    ! Finish

    return

  contains

    function dlog_ (y_a, y_b) result (dlog)

      real(WP), intent(in) :: y_a
      real(WP), intent(in) :: y_b
      real(WP)             :: dlog

      ! Calculate the logarithmic change between y_a and y_b

      if ((y_a > 0._WP .AND. y_b > 0._WP) .OR. &
          (y_a < 0._WP .AND. y_b < 0._WP)) then
         dlog = ABS(LOG10(y_b/y_a))
      else
         dlog = 0._WP
      endif

      ! Finish

      return

    end function dlog_

  end function dn_struct_
  
end module gyre_grid_factory
