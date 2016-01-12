! Module   : gyre_grid
! Purpose  : grid construction
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

module gyre_grid

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_grid_par
  use gyre_grid_weights
  use gyre_grid_util
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_grid

  ! Procedures

contains

  subroutine build_grid (ml, omega, gr_p, md_p, os_p, s, x, verbose)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, allocatable, intent(out)   :: s(:)
    real(WP), allocatable, intent(out)  :: x(:)
    logical, optional, intent(in)       :: verbose

    logical               :: write_info
    integer               :: n_s
    integer               :: s_
    real(WP), allocatable :: x_(:)

    if (PRESENT(verbose)) then
       write_info = verbose .AND. check_log_level('INFO')
    else
       write_info = check_log_level('DEBUG')
    endif

    ! Build a grid

    if (write_info) then

       write(OUTPUT_UNIT, 100) 'Building x grid'
100    format(A)

    endif

    ! Create the base grid, segment-by-segment

    n_s = ml%n_s

    allocate(x(0))
    allocate(s(0))

    seg_loop : do s_ = 1, ml%n_s

       x_ = seg_grid_(ml, s_, omega, gr_p, md_p, os_p)

       x = [x,x_]
       s = [s,SPREAD(s_, 1, SIZE(x_))]

       if (write_info) then
          write(OUTPUT_UNIT, 110) 'segment', s_, ':', SIZE(x_), 'points, x range', MINVAL(x_), '->', MAXVAL(x_)
110       format(3X,A,1X,I0,1X,A,1X,I0,1X,A,1X,F6.4,1X,A,1X,F6.4)
       endif

    end do seg_loop

    ! Add points at the center

    ! k_turn = 0
    ! x_turn = HUGE(0._WP)

    ! omega_loop : do j = 1, SIZE(omega)
    !    call find_turn(ml, s, x, omega(j), md_p, os_p, k_turn_j, x_turn_j)
    !    if (x_turn_j < x_turn) then
    !       k_turn = k_turn_j
    !       x_turn = x_turn_j
    !    endif
    ! end do omega_loop

    ! Finish

    return

  end subroutine build_grid

  !****

  function seg_grid_ (ml, s, omega, gr_p, md_p, os_p) result (x)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    real(WP), allocatable               :: x(:)

    real(WP), allocatable :: w(:)
    integer, allocatable  :: dn(:)
    integer               :: n_k
    real(WP), allocatable :: x_res(:)
    integer               :: i
    integer               :: j
    integer               :: k

    ! Construct the grid for the s'th segment

    ! First, make the base grid

    select case (gr_p%base_type)
    case ('MODEL')
       x = ml%x_base(s)
    case ('UNIFORM')
       w = uni_weights(gr_p%n_base)
       x = (1._WP-w)*ml%x_i(s) + w*ml%x_o(s)
    case ('GEOM')
       w = geo_weights(gr_p%n_base, gr_p%delta_base)
       x = (1._WP-w)*ml%x_i(s) + w*ml%x_o(s)
    case ('LOG')
       w = log_weights(gr_p%n_base, gr_p%delta_base)
       x = (1._WP-w)*ml%x_i(s) + w*ml%x_o(s)
    case default
       $ABORT(Invalid base_type)
    end select

    ! Now resample it by adding dn points to each subinterval

    ! First, determine dn

    dn = MAX(dn_dispersion_(ml, s, x, omega, gr_p, md_p, os_p), &
             dn_thermal_(ml, s, x, omega, gr_p, md_p, os_p), &
             dn_struct_(ml, s, x, gr_p))

    ! Add points

    n_k = SIZE(x)

    allocate(x_res(n_k + SUM(dn)))

    j = 1

    do k = 1, n_k-1
       do i = 1, dn(k)+1
          x_res(j) = x(k) + (i-1)*(x(k+1)-x(k))/(dn(k)+1)
          j = j + 1
       end do
    end do
    
    x_res(j) = x(n_k)

    call MOVE_ALLOC(x_res, x)

    ! Finish

    return

  end function seg_grid_

  !****

  function dn_dispersion_ (ml, s, x, omega, gr_p, md_p, os_p) result (dn)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s
    real(WP), intent(in)                :: x(:)
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer                             :: dn(SIZE(x)-1)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_k
    real(WP)                    :: beta_r_max(SIZE(x))
    real(WP)                    :: beta_i_max(SIZE(x))
    integer                     :: k
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
    real(WP)                    :: dphi_osc
    real(WP)                    :: dphi_exp

    ! Determine how many points dn to add to each cell of the grid x,
    ! such there are at least alpha_osc points per oscillatory
    ! wavelength and alpha_exp points per exponential wavelength
    !
    ! Wavelengths are calculated based on a local dispersion analysis
    ! of the adibatic/Cowling wave equation, for inertial frequencies
    ! specified by omega

    if (gr_p%alpha_osc > 0._WP .OR. gr_p%alpha_exp > 0._WP) then

       allocate(rt, SOURCE=r_rot_t(ml, md_p, os_p))

       ! At each point, determine the maximum absolute value of the
       ! real and imaginary parts of the local radial wavenumber beta,
       ! for all possible omega values

       n_k = SIZE(x)

       beta_r_max(1) = 0._WP
       beta_i_max(1) = 0._WP

       wavenumber_loop : do k = 2, n_k-1

          V_g = ml%V_2(s, x(k))*x(k)**2/ml%Gamma_1(s, x(k))
          As = ml%As(s, x(k))
          U = ml%U(s, x(k))
          c_1 = ml%c_1(s, x(k))

          beta_r_max(k) = 0._WP
          beta_i_max(k) = 0._WP

          omega_loop : do j = 1, SIZE(omega)

             omega_c = rt%omega_c(s, x(k), omega(j))

             lambda = rt%lambda(s, x(k), omega(j))
             l_i = rt%l_i(omega(j))
            
             ! Calculate the propagation discriminant gamma

             g_4 = -4._WP*V_g*c_1
             g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
             g_0 = -4._WP*lambda*As/c_1

             gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

             ! Update the wavenumber maxima

             if (gamma < 0._WP) then
               
                ! Propagation zone

                beta_r_max(k) = MAX(beta_r_max(k), ABS(0.5_WP*SQRT(-gamma))/x(k))
                beta_i_max(k) = MAX(beta_i_max(k), ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i))/x(k))

             else

                ! Evanescent zone

                beta_i_max(k) = MAX(beta_i_max(k), &
                     ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i - SQRT(gamma)))/x(k), &
                     ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i + SQRT(gamma)))/x(k))
             
             end if

          end do omega_loop

       end do wavenumber_loop

       beta_r_max(n_k) = 0._WP
       beta_i_max(n_k) = 0._WP

       ! Set up dn

       cell_loop : do k = 1, n_k-1

          ! Calculate the oscillatory and exponential phase change across
          ! the cell

          dphi_osc = MAX(beta_r_max(k), beta_r_max(k+1))*(x(k+1) - x(k))
          dphi_exp = MAX(beta_i_max(k), beta_i_max(k+1))*(x(k+1) - x(k))

          ! Set up dn

          dn(k) = MAX(FLOOR((gr_p%alpha_osc*dphi_osc)/TWOPI), FLOOR((gr_p%alpha_exp*dphi_exp)/TWOPI))

       end do cell_loop

    else

       dn = 0

    endif

    ! Finish

    return

  end function dn_dispersion_

  !****

  function dn_thermal_ (ml, s, x, omega, gr_p, md_p, os_p) result (dn)

    class(model_t), pointer, intent(in)  :: ml
    integer, intent(in)                  :: s
    real(WP), intent(in)                 :: x(:)
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer                              :: dn(SIZE(x)-1)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_k
    real(WP)                    :: beta_t_max(SIZE(x))
    integer                     :: k
    real(WP)                    :: V
    real(WP)                    :: nabla
    real(WP)                    :: c_rad
    real(WP)                    :: c_thm
    integer                     :: j
    real(WP)                    :: omega_c
    real(WP)                    :: dphi_thm

    ! Determine how many points dn to add to each cell of the grid x,
    ! such there are at least alpha_thm points per thermal length

    if (gr_p%alpha_thm > 0._WP) then

       allocate(rt, SOURCE=r_rot_t(ml, md_p, os_p))

       ! At each point, determine the maximum absolute value of the
       ! local thermal wavenumber beta_t, for all possible omega
       ! values

       n_k = SIZE(x)

       beta_t_max(1) = 0._WP

       wavenumber_loop : do k = 2, n_k-1

          V = ml%V_2(s, x(k))*x(k)**2
          nabla = ml%nabla(s, x(k))

          c_rad = ml%c_rad(s, x(k))
          c_thm = ml%c_thm(s, x(k))

          beta_t_max(k) = 0._WP

          omega_loop : do j = 1, SIZE(omega)

             omega_c = rt%omega_c(s, x(k), omega(j))

             beta_t_max(k) = MAX(beta_t_max(k), SQRT(ABS(V*nabla*omega_c*c_thm/c_rad))/x(k))

          end do omega_loop

       end do wavenumber_loop

       beta_t_max(n_k) = 0._WP

       ! Set up dn
       
       cell_loop : do k = 1, n_k-1

          ! Calculate the thermal phase change across the cell

          dphi_thm = MAX(beta_t_max(k), beta_t_max(k+1))*(x(k+1) - x(k))

          ! Set up dn

          dn(k) = FLOOR((gr_p%alpha_thm*dphi_thm)/TWOPI)

       end do cell_loop

    else

       dn = 0

    endif

    ! Finish

    return

    return

  end function dn_thermal_

  !****

  function dn_struct_ (ml, s, x, gr_p) result (dn)

    class(model_t), pointer, intent(in)  :: ml
    integer, intent(in)                  :: s
    real(WP), intent(in)                 :: x(:)
    type(grid_par_t), intent(in)         :: gr_p
    integer                              :: dn(SIZE(x)-1)

    integer :: k

    ! Determine how many points dn to add to each cell of the grid x,
    ! such there are at least alpha_str points per dex change in the
    ! structure variables (V, As, Gamma_1, c_1, & U)

    if (gr_p%alpha_str > 0) then

       ! Set up dn

       cell_loop : do k = 1, SIZE(x)-1

          dn(k) = FLOOR(gr_p%alpha_str*dlog_(ml%V_2(s, x(k)), ml%V_2(s, x(k+1)))) + &
                  FLOOR(gr_p%alpha_str*dlog_(ml%As(s, x(k)), ml%As(s, x(k+1)))) + &
                  FLOOR(gr_p%alpha_str*dlog_(ml%Gamma_1(s, x(k)), ml%Gamma_1(s, x(k+1)))) + &
                  FLOOR(gr_p%alpha_str*dlog_(ml%c_1(s, x(k)), ml%c_1(s, x(k+1)))) + &
                  FLOOR(gr_p%alpha_str*dlog_(ml%U(s, x(k)), ml%U(s, x(k+1))))

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

  ! !****

  ! function dn_center_ (ml, s, x, omega, gr_p, md_p, os_p) result (dn)

  !   class(model_t), pointer,  intent(in) :: ml
  !   integer, intent(in)                  :: s
  !   real(WP), intent(in)                 :: x(:)
  !   real(WP), intent(in)                 :: omega(:)
  !   type(grid_par_t), intent(in)         :: gr_p
  !   type(mode_par_t), intent(in)         :: md_p
  !   type(osc_par_t), intent(in)          :: os_p
  !   integer                              :: dn(SIZE(x)-1)
    
  !   integer  :: k_turn
  !   real(WP) :: x_turn
  !   integer  :: j
  !   integer  :: k_turn_j
  !   real(WP) :: x_turn_j
  !   integer  :: n_add
  !   integer  :: n_k
  !   integer  :: k

  !   ! Determine how many points dn to add to each cell of the grid x,
  !   ! such there are at least n_center points covering the evanescent
  !   ! region at the center
  !   !
  !   ! Evanescence is determined based on a local dispersion analysis
  !   ! of the adibatic/Cowling wave equation, for all possible omega
  !   ! values

  !   if (gr_p%n_center > 0) then

  !      ! Locate the innermost turning point

  !      k_turn = 0
  !      x_turn = HUGE(0._WP)

  !      omega_loop : do j = 1, SIZE(omega)
  !         call find_turn(ml, s, x, omega(j), md_p, os_p, k_turn_j, x_turn_j)
  !         if (x_turn_j < x_turn) then
  !            k_turn = k_turn_j
  !            x_turn = x_turn_j
  !         endif
  !      end do omega_loop

  !      ! Determine how many points need to be added

  !      n_add = MAX(gr_p%n_center-k_turn, 0)

  !      ! Set up dn

  !      cell_loop : do k = 1, n_k-1
  !         if (k <= k_turn) then
  !            dn(k) = CEILING(n_add*(x(k_turn+1)/x_turn)/k_turn)
  !         else
  !            dn(k) = 0
  !         endif
  !      end do cell_loop

  !   else

  !      dn = 0

  !   endif

  !   ! Finish

  !   return

  ! end function dn_center_

end module gyre_grid
