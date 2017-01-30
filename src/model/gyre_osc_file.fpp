! Module   : gyre_osc_file
! Purpose  : read OSC files
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

module gyre_osc_file

  ! Uses

  use core_kinds
  use core_order

  use gyre_constants
  use gyre_evol_model
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_osc_model

  ! Procedures

contains

  subroutine read_osc_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    character(:), allocatable   :: data_format
    integer                     :: unit
    integer                     :: n
    integer                     :: iconst
    integer                     :: ivar
    integer                     :: iabund
    integer                     :: ivers
    real(WP), allocatable       :: global_data(:)
    real(WP), allocatable       :: point_data(:,:)
    integer                     :: i
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: p(:)
    real(WP), allocatable       :: rho(:) 
    real(WP), allocatable       :: T(:) 
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kap(:)
    real(WP), allocatable       :: kap_rho(:)
    real(WP), allocatable       :: kap_T(:)
    real(WP), allocatable       :: eps(:)
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thm(:)
    real(WP), allocatable       :: c_dif(:)
    real(WP), allocatable       :: c_eps_ad(:)
    real(WP), allocatable       :: c_eps_S(:)
    real(WP), allocatable       :: kap_ad(:)
    real(WP), allocatable       :: kap_S(:)
    type(evol_model_t), pointer :: em

    ! Open the OSC-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from OSC file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif
          
    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, iabund, ivers

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'File version', ivers
120    format(3X,A,1X,F4.2,1X,A)
    endif

    ! Read the data

    if (ml_p%data_format /= '') then
       data_format = ml_p%data_format
    else
       data_format = '(1P5E19.12)'
    endif

    allocate(global_data(iconst))
    allocate(point_data(ivar+iabund,n))

    read(unit, data_format) global_data

    read_loop : do i = 1,n
       read(unit, data_format) point_data(:,i)
    end do read_loop

    close(unit)

    point_data = point_data(:,n:1:-1)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', n, 'points'
130    format(3X,A,1X,I0,1X,A)
    endif
    
    ! Extract structure data

    M_star = global_data(1)
    R_star = global_data(2)
    L_star = global_data(3)

    r = point_data(1,:)
    m = EXP(point_data(2,:))

    T = point_data(3,:)
    P = point_data(4,:)
    rho = point_data(5,:)

    Gamma_1 = point_data(10,:)
    nabla_ad = point_data(11,:)
    delta = point_data(12,:)

    As = point_data(15,:)

    nabla = point_data(6,:)
    kap = point_data(8,:)
    kap_T = point_data(17,:)
    kap_rho = point_data(18,:)
    eps = point_data(9,:)
    eps_T = point_data(19,:)
    eps_rho = point_data(20,:)

    Omega_rot = point_data(16,:)

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, m)

    ! Calculate dimensionless structure data

    allocate(V_2(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(P*r*x**2)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

    kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
    kap_S = kap_T - delta*kap_rho

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kap*rho*L_star)
    c_thm = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_dif = (kap_ad-4._WP*nabla_ad)*V_2*x**2*nabla + V_2*x**2*nabla_ad

    c_eps_ad = 4._WP*PI*rho*(nabla_ad*eps_T + eps_rho/Gamma_1)*R_star**3/L_star
    c_eps_S = 4._WP*PI*rho*(eps_T - delta*eps_rho)*R_star**3/L_star

    if (ml_p%uniform_rot) then
       allocate(Omega_rot(n))
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%define(I_V_2, V_2)
    call em%define(I_AS, As)
    call em%define(I_U, U)
    call em%define(I_C_1, c_1)

    call em%define(I_GAMMA_1, Gamma_1)
    call em%define(I_DELTA, delta)
    call em%define(I_NABLA_AD, nabla_ad)
    call em%define(I_NABLA, nabla)
    call em%define(I_BETA_RAD, beta_rad)

    call em%define(I_C_RAD, c_rad)
    call em%define(I_C_THM, c_thm)
    call em%define(I_C_DIF, c_dif)
    call em%define(I_C_EPS_AD, c_eps_ad)
    call em%define(I_C_EPS_S, c_eps_S)
    call em%define(I_KAP_AD, kap_ad)
    call em%define(I_KAP_S, kap_S)

    call em%define(I_OMEGA_ROT, Omega_rot)

    ! Return a pointer

    ml => em

    ! Finish

    return

  end subroutine read_osc_model

end module gyre_osc_file
