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

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_evol_model
  use gyre_ext
  use gyre_input
  use gyre_grid
  use gyre_grid_par
  use gyre_mesa_file
  use gyre_mode
  use gyre_model
  use gyre_mode_par
  use gyre_nad_bvp
  use gyre_osc_par
  use gyre_num_par
  use gyre_rad_bvp
  use gyre_scan_par
  use gyre_search
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

    if (ASSOCIATED(ml_m)) deallocate(ml_m)
    if (ALLOCATED(x_ml_m)) deallocate(x_ml_m)

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

    if (ASSOCIATED(ml_m)) deallocate(ml_m)

    call read_mesa_model(file, deriv_type, .TRUE., ec, x_ml_m)

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
                                       deriv_type, add_center=add_center))

    if(add_center) then
       x_ml_m = [0._WP,r/R_star]
    else
       x_ml_m = r/R_star
    endif

    ! Finish

    return

  end subroutine gyre_set_model

!****

  subroutine gyre_get_modes (l, filename, non_ad, user_sub, ipar, rpar)

    integer, intent(in)          :: l
    character(LEN=*), intent(in) :: filename
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

    integer                       :: unit
    type(mode_par_t), allocatable :: mp(:)
    type(osc_par_t), allocatable  :: op(:)
    type(num_par_t), allocatable  :: np(:)
    type(scan_par_t), allocatable :: sp(:)
    type(grid_par_t), allocatable :: shoot_gp(:)
    type(grid_par_t), allocatable :: recon_gp(:)
    integer                       :: n_md
    integer                       :: d_md
    type(mode_t), allocatable     :: md(:)
    integer                       :: i
    type(osc_par_t), allocatable  :: op_sel(:)
    type(num_par_t), allocatable  :: np_sel(:)
    type(grid_par_t), allocatable :: shoot_gp_sel(:)
    type(grid_par_t), allocatable :: recon_gp_sel(:)
    type(scan_par_t), allocatable :: sp_sel(:)
    real(WP)                      :: x_i
    real(WP)                      :: x_o
    real(WP), allocatable         :: omega(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max
    real(WP), allocatable         :: x_sh(:)
    class(r_bvp_t), allocatable   :: ad_bp
    class(c_bvp_t), allocatable   :: nad_bp

    $ASSERT(ASSOCIATED(ml_m),No model provided)

    ! Read parameters

    open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

    call read_mode_par(unit, mp)
    call read_osc_par(unit, op)
    call read_num_par(unit, np)
    call read_shoot_grid_par(unit, shoot_gp)
    call read_recon_grid_par(unit, recon_gp)
    call read_scan_par(unit, sp)

    close(unit)

    ! Loop through modepars

    d_md = 128
    n_md = 0

    allocate(md(d_md))

    mp_loop : do i = 1, SIZE(mp)

       if (mp(i)%l == l) then

          ! Select parameters according to tags

          call select_par(op, mp(i)%tag, op_sel, last=.TRUE.)
          call select_par(np, mp(i)%tag, np_sel, last=.TRUE.)
          call select_par(shoot_gp, mp(i)%tag, shoot_gp_sel)
          call select_par(recon_gp, mp(i)%tag, recon_gp_sel)
          call select_par(sp, mp(i)%tag, sp_sel)

          $ASSERT(SIZE(op_sel) == 1,No matching osc parameters)
          $ASSERT(SIZE(np_sel) == 1,No matching num parameters)
          $ASSERT(SIZE(shoot_gp_sel) >= 1,No matching shoot_grid parameters)
          $ASSERT(SIZE(recon_gp_sel) >= 1,No matching recon_grid parameters)
          $ASSERT(SIZE(sp_sel) >= 1,No matching scan parameters)

          ! Set up the frequency array

          x_i = x_ml_m(1)
          x_o = x_ml_m(SIZE(x_ml_m))

          call build_scan(sp_sel, ml_m, mp(i), op_sel(1), x_i, x_o, omega)

          ! Set up the shooting grid
          
          call build_grid(shoot_gp_sel, ml_m, mp(i), op_sel(1), omega, x_ml_m, x_sh)

          if (np_sel(1)%restrict_roots) then
             omega_min = MINVAL(omega)
             omega_max = MAXVAL(omega)
          else
             omega_min = -HUGE(0._WP)
             omega_max = HUGE(0._WP)
          endif

          ! Set up the bvp's

          if(mp(i)%l == 0 .AND. op_sel(1)%reduce_order) then
             allocate(ad_bp, SOURCE=rad_bvp_t(x_sh, ml_m, mp(i), op_sel(1), np_sel(1), omega_min, omega_max))
          else
             allocate(ad_bp, SOURCE=ad_bvp_t(x_sh, ml_m, mp(i), op_sel(1), np_sel(1), omega_min, omega_max))
          endif

          if (non_ad) then
             allocate(nad_bp, SOURCE=nad_bvp_t(x_sh, ml_m, mp(i), op_sel(1), np_sel(1), omega_min, omega_max))
          endif

          ! Find modes

          if (non_ad) then
             n_md = 0
             call scan_search(ad_bp, np_sel(1), omega, process_root_ad)
             call prox_search(nad_bp, mp(i), np_sel(1), op_sel(1), md(:n_md), process_root_nad)
          else
             call scan_search(ad_bp, np_sel(1), omega, process_root_ad)
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

    subroutine process_root_ad (omega, n_iter, discrim_ref)

      real(WP), intent(in)      :: omega
      integer, intent(in)       :: n_iter
      type(r_ext_t), intent(in) :: discrim_ref

      real(WP), allocatable :: x_rc(:)
      integer               :: n
      real(WP)              :: x_ref
      real(WP), allocatable :: y(:,:)
      real(WP)              :: y_ref(6)
      type(r_ext_t)         :: discrim
      type(mode_t)          :: md_new
      integer               :: retcode

      ! Build the reconstruction grid

      call build_grid(recon_gp_sel, ml_m, mp(i), op_sel(1), [omega], x_sh, x_rc, verbose=.FALSE.)

      ! Reconstruct the solution

      x_ref = MIN(MAX(op_sel(1)%x_ref, x_sh(1)), x_sh(SIZE(x_sh)))

      n = SIZE(x_rc)

      allocate(y(6,n))

      call ad_bp%recon(omega, x_rc, x_ref, y, y_ref, discrim)

      ! Create the mode

      md_new = mode_t(ml_m, mp(i), op_sel(1), CMPLX(omega, KIND=WP), c_ext_t(discrim), &
                      x_rc, CMPLX(y, KIND=WP), x_ref, CMPLX(y_ref, KIND=WP))

      if (md_new%n_pg < mp(i)%n_pg_min .OR. md_new%n_pg > mp(i)%n_pg_max) return

      md_new%n_iter = n_iter
      md_new%chi = ABS(discrim)/ABS(discrim_ref)

      ! Store it or process it

      if (non_ad) then

         n_md = n_md + 1

         if (n_md > d_md) then
            d_md = 2*d_md
            call reallocate(md, [d_md])
         endif
         
         md(n_md) = md_new

         call md(n_md)%prune()

      else

         retcode = 0

         call user_sub(md_new, ipar, rpar, retcode)

      endif

      ! Finish

      return

    end subroutine process_root_ad

!****

    subroutine process_root_nad (omega, n_iter, discrim_ref)

      complex(WP), intent(in)   :: omega
      integer, intent(in)       :: n_iter
      type(r_ext_t), intent(in) :: discrim_ref

      real(WP), allocatable    :: x_rc(:)
      integer                  :: n
      real(WP)                 :: x_ref
      complex(WP), allocatable :: y(:,:)
      complex(WP)              :: y_ref(6)
      type(c_ext_t)            :: discrim
      type(mode_t)             :: md_new
      integer                  :: retcode

      ! Build the reconstruction grid

      call build_grid(recon_gp_sel, ml_m, mp(i), op_sel(1), [REAL(omega)], x_sh, x_rc, verbose=.FALSE.)

      ! Reconstruct the solution

      x_ref = MIN(MAX(op_sel(1)%x_ref, x_sh(1)), x_sh(SIZE(x_sh)))

      n = SIZE(x_rc)

      allocate(y(6,n))

      call nad_bp%recon(omega, x_rc, x_ref, y, y_ref, discrim)

      ! Create the mode

      md_new = mode_t(ml_m, mp(i), op_sel(1), CMPLX(omega, KIND=WP), discrim, &
           x_rc, y, x_ref, y_ref)

      md_new%n_iter = n_iter
      md_new%chi = ABS(discrim)/ABS(discrim_ref)

      if (check_log_level('INFO')) then
         write(OUTPUT_UNIT, 120) md_new%mp%l, md_new%n_pg, md_new%n_p, md_new%n_g, &
              md_new%omega, real(md_new%chi), md_new%n_iter, md_new%n
120      format(4(2X,I8),3(2X,E24.16),2X,I6,2X,I7)
      endif

      ! Process it

      retcode = 0

      call user_sub(md_new, ipar, rpar, retcode)

      ! Finish

      return

    end subroutine process_root_nad

  end subroutine gyre_get_modes

end module gyre_lib
