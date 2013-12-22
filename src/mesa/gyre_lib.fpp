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

  use gyre_bvp
  use gyre_bvp_ad
  use gyre_bvp_rad
  use gyre_coeffs
  use gyre_coeffs_evol
  use gyre_mesa_file
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

  type(coeffs_evol_t), allocatable, save :: cf_m
  real(WP), allocatable, save            :: x_cf_m(:)

  ! Access specifiers

  private

  public :: WP
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

    if(ALLOCATED(cf_m)) then
       call cf_m%final()
       deallocate(cf_m)
    endif

    if(ALLOCATED(x_cf_m)) deallocate(x_cf_m)

    call final_parallel()

    ! Finish

    return

  end subroutine gyre_final

!****

  subroutine gyre_read_model (file, G, deriv_type)

    character(LEN=*), intent(in) :: file
    real(WP), intent(in)         :: G
    character(LEN=*), intent(in) :: deriv_type

    ! Allocate the model

    if(ALLOCATED(cf_m)) then
       call cf_m%final()
       deallocate(cf_m)
    endif

    allocate(coeffs_evol_t::cf_m)

    ! Read the model

    call read_mesa_file(file, G, deriv_type, cf_m, x_cf_m)

    ! Finish

    return

  end subroutine gyre_read_model
  
!****

  subroutine gyre_set_model (G, M_star, R_star, L_star, r, w, p, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, &
                             Omega_rot, deriv_type)

    real(WP), intent(in)         :: G
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

    if(ALLOCATED(cf_m)) then
       call cf_m%final()
       deallocate(cf_m)
    endif

    allocate(coeffs_evol_t::cf_m)

    ! Set the model by storing coefficients

    m = w/(1._WP+w)*M_star

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    call cf_m%init(G, M_star, R_star, L_star, r, m, p, rho, T, N2, &
                   Gamma_1, nabla_ad, delta, Omega_rot, &
                   nabla, kappa, kappa_rho, kappa_T, &
                   epsilon, epsilon_rho, epsilon_T, &
                   deriv_type, add_center)

    if(add_center) then
       x_cf_m = [0._WP,r/R_star]
    else
       x_cf_m = r/R_star
    endif

    ! Finish

    return

  end subroutine gyre_set_model

!****

  subroutine gyre_get_modes (file, user_sub, ipar, rpar, dummy)

    character(LEN=*), intent(in) :: file
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
    logical, intent(in), optional :: dummy

    integer                      :: unit
    type(oscpar_t), allocatable  :: op(:)
    type(numpar_t), allocatable  :: np(:)
    type(scanpar_t), allocatable :: sp(:)
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    integer                      :: i
    type(numpar_t), allocatable  :: np_sel(:)
    type(gridpar_t), allocatable :: shoot_gp_sel(:)
    type(gridpar_t), allocatable :: recon_gp_sel(:)
    type(scanpar_t), allocatable :: sp_sel(:)
    real(WP), allocatable        :: omega(:)
    class(bvp_t), allocatable    :: bp
    type(mode_t), allocatable    :: md(:)
    integer                      :: j
    integer                      :: retcode

    $ASSERT(ALLOCATED(cf_m),No model provided)

    ! Read parameters

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_oscpar(unit, op)
    call read_numpar(unit, np)
    call read_shoot_gridpar(unit, shoot_gp)
    call read_recon_gridpar(unit, recon_gp)
    call read_scanpar(unit, sp)

    close(unit)

    ! Loop through oscpars

    op_loop : do i = 1, SIZE(op)

       ! Select parameters according to tags

       call select_par(np, op(i)%tag, np_sel, last=.TRUE.)
       call select_par(shoot_gp, op(i)%tag, shoot_gp_sel)
       call select_par(recon_gp, op(i)%tag, recon_gp_sel)
       call select_par(sp, op(i)%tag, sp_sel)

       $ASSERT(SIZE(np_sel) == 1,No matching num parameters)
       $ASSERT(SIZE(shoot_gp_sel) >= 1,No matching shoot_grid parameters)
       $ASSERT(SIZE(recon_gp_sel) >= 1,No matching recon_grid parameters)
       $ASSERT(SIZE(sp_sel) >= 1,No matching scan parameters)

       ! Set up the frequency array

       call build_scan(sp_sel, cf_m, op(i), shoot_gp_sel, x_cf_m, omega)

       ! Store the frequency range in shoot_gp_sel

       shoot_gp_sel%omega_a = MINVAL(omega)
       shoot_gp_sel%omega_b = MAXVAL(omega)

       ! Set up bp

       if(op(i)%l == 0 .AND. np_sel(1)%reduce_order) then
          allocate(bvp_rad_t::bp)
       else
          allocate(bvp_ad_t::bp)
       endif

       call bp%init(cf_m, op(i), np_sel(1), shoot_gp_sel, recon_gp_sel, x_cf_m)

       ! Find modes

       call scan_search(bp, omega, md)

       $if($GFORTRAN_PR57922)
       select type (bp)
       type is (bvp_rad_t)
          call bp%final()
       type is (bvp_ad_t)
          call bp%final()
       class default
          $ABORT(Invalid type)
       end select
       $endif

       if(ALLOCATED(bp)) deallocate(bp)

       ! Process the modes

       retcode = 0

       mode_loop : do j = 1,SIZE(md)
          if(retcode == 0) then
             call user_sub(md(j), ipar, rpar, retcode)
          endif
          $if($GFORTRAN_PR57922)
          call md(j)%final()
          $endif
       end do mode_loop

       ! Loop around

    end do op_loop

    ! Finish

    return

  end subroutine gyre_get_modes

end module gyre_lib
