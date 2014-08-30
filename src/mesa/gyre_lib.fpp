! Module   : gyre_lib
! Purpose  : library interface for use in MESA
!
! Copyright 2013 Rich Townsend
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

module gyre_lib

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants
  use gyre_bvp
  use gyre_ad_bvp
  use gyre_rad_bvp
  use gyre_nad_bvp
  use gyre_model
  use gyre_evol_model
  use gyre_mesa_file
  use gyre_modepar
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_scanpar
  use gyre_search
  use gyre_mode
  use gyre_input
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  class(model_t), pointer, save :: ml_m => null()
  real(WP), allocatable, save   :: x_ml_m(:)

  ! Access specifiers

  private

  public :: WP

  public :: G_GRAVITY
  public :: C_LIGHT
  public :: A_RADIATION
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN

  public :: mode_t
  public :: gyre_init
  public :: gyre_read_model
  public :: gyre_set_model
  public :: gyre_get_modes

  ! Procedures

contains

  subroutine gyre_init ()

    ! Initialize

    call init_parallel()

    call set_log_level('WARN')

    ! Finish

    return

  end subroutine gyre_init

!****

  subroutine gyre_final()

    ! Finalize

    if(ASSOCIATED(ml_m)) deallocate(ml_m)
    if(ALLOCATED(x_ml_m)) deallocate(x_ml_m)

    call final_parallel()

    ! Finish

    return

  end subroutine gyre_final

!****

  subroutine gyre_read_model (file, deriv_type)

    character(LEN=*), intent(in) :: file
    character(LEN=*), intent(in) :: deriv_type

    type(evol_model_t) :: ec

    ! Read the model

    if(ASSOCIATED(ml_m)) deallocate(ml_m)

    call read_mesa_model(file, deriv_type, .FALSE., ec, x_ml_m)

    allocate(ml_m, SOURCE=ec)

    ! Finish

    return

  end subroutine gyre_read_model
  
!****

  subroutine gyre_set_model (M_star, R_star, L_star, r, w, p, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, &
                             Omega_rot, deriv_type)

    real(WP), intent(in)         :: M_star
    real(WP), intent(in)         :: R_star
    real(WP), intent(in)         :: L_star
    real(WP), intent(in)         :: r(:)
    real(WP), intent(in)         :: w(:)
    real(WP), intent(in)         :: p(:)
    real(WP), intent(in)         :: rho(:)
    real(WP), intent(in)         :: T(:)
    real(WP), intent(in)         :: N2(:)
    real(WP), intent(in)         :: Gamma_1(:)
    real(WP), intent(in)         :: nabla_ad(:)
    real(WP), intent(in)         :: delta(:)
    real(WP), intent(in)         :: nabla(:)
    real(WP), intent(in)         :: kappa(:)
    real(WP), intent(in)         :: kappa_rho(:)
    real(WP), intent(in)         :: kappa_T(:)
    real(WP), intent(in)         :: epsilon(:)
    real(WP), intent(in)         :: epsilon_rho(:)
    real(WP), intent(in)         :: epsilon_T(:)
    real(WP), intent(in)         :: Omega_rot(:)
    character(LEN=*), intent(in) :: deriv_type

    real(WP), allocatable :: m(:)
    logical               :: add_center

    ! Allocate the model

    if(ASSOCIATED(ml_m)) deallocate(ml_m)
    allocate(evol_model_t::ml_m)

    ! Set the model by storing coefficients

    m = w/(1._WP+w)*M_star

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    allocate(ml_m, SOURCE=evol_model_t(M_star, R_star, L_star, r, m, p, rho, T, N2, &
                                       Gamma_1, nabla_ad, delta, Omega_rot, &
                                       nabla, kappa, kappa_rho, kappa_T, &
                                       epsilon, epsilon_rho, epsilon_T, &
                                       deriv_type, add_center))

    if(add_center) then
       x_ml_m = [0._WP,r/R_star]
    else
       x_ml_m = r/R_star
    endif

    ! Finish

    return

  end subroutine gyre_set_model

!****

  subroutine gyre_get_modes (l, file, non_ad, user_sub, ipar, rpar)

    integer, intent(in)          :: l
    character(LEN=*), intent(in) :: file
    logical, intent(in)          :: non_ad
    interface
       subroutine user_sub (md, ipar, rpar, retcode)
         import mode_t
         import WP
         type(mode_t), intent(in) :: md
         integer, intent(inout)   :: ipar(:)
         real(WP), intent(inout)  :: rpar(:)
         integer, intent(out)     :: retcode
       end subroutine user_sub
    end interface
    integer, intent(inout)  :: ipar(:)
    real(WP), intent(inout) :: rpar(:)

    integer                      :: unit
    type(modepar_t), allocatable :: mp(:)
    type(oscpar_t), allocatable  :: op(:)
    type(numpar_t), allocatable  :: np(:)
    type(scanpar_t), allocatable :: sp(:)
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    integer                      :: n_md
    integer                      :: d_md
    type(mode_t), allocatable    :: md(:)
    integer                      :: i
    type(oscpar_t), allocatable  :: op_sel(:)
    type(numpar_t), allocatable  :: np_sel(:)
    type(gridpar_t), allocatable :: shoot_gp_sel(:)
    type(gridpar_t), allocatable :: recon_gp_sel(:)
    type(scanpar_t), allocatable :: sp_sel(:)
    integer                      :: n_op_sel
    integer                      :: n_np_sel
    real(WP), allocatable        :: omega(:)
    class(bvp_t), allocatable    :: ad_bp
    class(bvp_t), allocatable    :: nad_bp

    $ASSERT(ASSOCIATED(ml_m),No model provided)

    ! Read parameters

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_modepar(unit, mp)
    call read_oscpar(unit, op)
    call read_numpar(unit, np)
    call read_shoot_gridpar(unit, shoot_gp)
    call read_recon_gridpar(unit, recon_gp)
    call read_scanpar(unit, sp)

    close(unit)

    ! Loop through modepars

    d_md = 128
    n_md = 0

    allocate(md(d_md))

    mp_loop : do i = 1, SIZE(mp)

       if (mp(i)%l == l) then

          ! Select parameters according to tags

          call select_par(op, mp(i)%tag, op_sel)
          call select_par(np, mp(i)%tag, np_sel)
          call select_par(shoot_gp, mp(i)%tag, shoot_gp_sel)
          call select_par(recon_gp, mp(i)%tag, recon_gp_sel)
          call select_par(sp, mp(i)%tag, sp_sel)

          n_op_sel = SIZE(op_sel)
          n_np_sel = SIZE(np_sel)

          $ASSERT(n_op_sel >= 1,No matching osc parameters)
          $ASSERT(n_np_sel >= 1,No matching num parameters)
          $ASSERT(SIZE(shoot_gp_sel) >= 1,No matching shoot_grid parameters)
          $ASSERT(SIZE(recon_gp_sel) >= 1,No matching recon_grid parameters)
          $ASSERT(SIZE(sp_sel) >= 1,No matching scan parameters)

          ! Set up the frequency array

          call build_scan(sp_sel, ml_m, mp(i), op_sel(1), shoot_gp_sel, x_ml_m, omega)

          ! Store the frequency range in shoot_gp_sel

          shoot_gp_sel%omega_a = MINVAL(omega)
          shoot_gp_sel%omega_b = MAXVAL(omega)

          ! Set up the bvp's

          if(mp(i)%l == 0 .AND. op_sel(1)%reduce_order) then
             allocate(ad_bp, SOURCE=rad_bvp_t(ml_m, mp(i), op_sel(1), np_sel(1), shoot_gp_sel, recon_gp_sel, x_ml_m))
          else
             allocate(ad_bp, SOURCE=ad_bvp_t(ml_m, mp(i), op_sel(1), np_sel(1), shoot_gp_sel, recon_gp_sel, x_ml_m))
          endif

          if (non_ad) then
             allocate(nad_bp, SOURCE=nad_bvp_t(ml_m, mp(i), op_sel(1), np_sel(1), shoot_gp_sel, recon_gp_sel, x_ml_m))
          endif

          ! Find modes

          if (non_ad) then
             n_md = 0
             call scan_search(ad_bp, np_sel(n_np_sel), omega, store_mode)
             call prox_search(nad_bp, np_sel(n_np_sel), md(:n_md), process_mode)
          else
             call scan_search(ad_bp, np_sel(n_np_sel), omega, process_mode)
          endif

          ! Clean up

          deallocate(ad_bp)
          if (non_ad) deallocate(nad_bp)

       end if

       ! Loop around

    end do mp_loop

    ! Finish

    return

  contains

    subroutine store_mode (md_new)

      type(mode_t), intent(in) :: md_new

      ! Store the mode

      n_md = n_md + 1

      if (n_md > d_md) then
         d_md = 2*d_md
         call reallocate(md, [d_md])
      endif

      md(n_md) = md_new

      call md(n_md)%prune()

      ! Finish

      return

    end subroutine store_mode

    subroutine process_mode (md_new)

      type(mode_t), intent(in) :: md_new

      integer :: retcode

      ! Process the mode

      retcode = 0

      call user_sub(md_new, ipar, rpar, retcode)

    end subroutine process_mode

  end subroutine gyre_get_modes

end module gyre_lib
