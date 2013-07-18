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
  use gyre_rad_bvp
  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_therm_coeffs
  use gyre_evol_therm_coeffs
  use gyre_mesa_file
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_ad_search
  use gyre_rad_search
  use gyre_mode
  use gyre_input
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  class(base_coeffs_t), allocatable, save  :: bc_m
  class(therm_coeffs_t), allocatable, save :: tc_m
  real(WP), allocatable, save              :: x_bc_m(:)

  ! Access specifiers

  private

  public :: WP
  public :: mode_t
  public :: gyre_init
  public :: gyre_read_model
  public :: gyre_set_model
  public :: gyre_get_modes
  public :: gyre_get_radial_modes

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

  subroutine gyre_read_model (file, G, deriv_type)

    character(LEN=*), intent(in) :: file
    real(WP), intent(in)         :: G
    character(LEN=*), intent(in) :: deriv_type

    ! Read the model

    call read_mesa_file(file, G, deriv_type, bc_m, tc_m, x_bc_m)

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

    ! Set the model by storing coefficients

    if(ALLOCATED(bc_m)) then
       $if($GFORTRAN_PR57922)
       call bc_m%final()
       $endif
       deallocate(bc_m)
    endif

    if(ALLOCATED(tc_m)) then
       $if($GFORTRAN_PR57922)
       call tc_m%final()
       $endif
       deallocate(tc_m)
    endif

    allocate(evol_base_coeffs_t::bc_m)
    allocate(evol_therm_coeffs_t::tc_m)

    m = w/(1._WP+w)*M_star

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    select type (bc_m)
    type is (evol_base_coeffs_t)
       call bc_m%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                      N2, Gamma_1, nabla_ad, delta, &
                      Omega_rot, deriv_type, add_center)
    class default
       $ABORT(Invalid bc_m type)
    end select

    select type (tc_m)
    type is (evol_therm_coeffs_t)
       call tc_m%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                      Gamma_1, nabla_ad, delta, nabla,  &
                      kappa, kappa_rho, kappa_T, &
                      epsilon, epsilon_rho, epsilon_T, deriv_type, add_center)
    class default
       $ABORT(Invalid tc_m type)
    end select

    if(add_center) then
       x_bc_m = [0._WP,r/R_star]
    else
       x_bc_m = r/R_star
    endif

    ! Finish

    return

  end subroutine gyre_set_model

!****

  subroutine gyre_get_modes (file, user_sub, ipar, rpar)

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

    integer                      :: unit
    type(oscpar_t)               :: op
    type(numpar_t)               :: np
    real(WP), allocatable        :: omega(:)
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    type(ad_bvp_t)               :: bp
    type(mode_t), allocatable    :: md(:)
    integer                      :: j
    integer                      :: retcode

    ! Read parameters

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_oscpar(unit, op)
    call read_numpar(unit, np)
    call read_shoot_gridpar(unit, shoot_gp)
    call read_recon_gridpar(unit, recon_gp)
    call read_scanpar(unit, bc_m, op, shoot_gp, x_bc_m, omega)

    close(unit)

    ! Set up bp

    call bp%init(bc_m, op, np, shoot_gp, recon_gp, x_bc_m, tc_m)

    ! Find modes

    call ad_scan_search(bp, omega, md)

    ! Process the modes

    mode_loop : do j = 1,SIZE(md)
       call user_sub(md(j), ipar, rpar, retcode)
       if(retcode /= 0) exit mode_loop
       $if($GFORTRAN_PR57922)
       call md(j)%final()
       $endif
    end do mode_loop

    ! Finish

    return

  end subroutine gyre_get_modes

!****

  subroutine gyre_get_radial_modes (file, user_sub, ipar, rpar)

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

    integer                      :: unit
    type(oscpar_t)               :: op
    type(numpar_t)               :: np
    real(WP), allocatable        :: omega(:)
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    type(rad_bvp_t)              :: bp
    type(mode_t), allocatable    :: md(:)
    integer                      :: j
    integer                      :: retcode

    ! Read parameters

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_oscpar(unit, op)
    call read_numpar(unit, np)
    call read_shoot_gridpar(unit, shoot_gp)
    call read_recon_gridpar(unit, recon_gp)
    call read_scanpar(unit, bc_m, op, shoot_gp, x_bc_m, omega)

    close(unit)

    ! Set up bp

    call bp%init(bc_m, op, np, shoot_gp, recon_gp, x_bc_m, tc_m)

    ! Find modes

    call rad_scan_search(bp, omega, md)

    ! Process the modes

    mode_loop : do j = 1,SIZE(md)
       call user_sub(md(j), ipar, rpar, retcode)
       if(retcode /= 0) exit mode_loop
    end do mode_loop

    ! Finish

    return

  end subroutine gyre_get_radial_modes

end module gyre_lib
