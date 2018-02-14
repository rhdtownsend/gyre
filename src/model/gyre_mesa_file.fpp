! Module   : gyre_mesa_file
! Purpose  : read MESA files
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

module gyre_mesa_file

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

  ! Parameters

  integer, parameter :: N_COLS_V0_01 = 18
  integer, parameter :: N_COLS_V0_19 = 18
  integer, parameter :: N_COLS_V1_0X = 18

  ! Access specifiers

  private

  public :: read_mesa_model
  public :: init_mesa_model

  ! Procedures

contains

    subroutine read_mesa_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                     :: unit
    integer                     :: n
    real(WP)                    :: global_data(3)
    integer                     :: version
    integer                     :: n_cols
    real(WP), allocatable       :: point_data(:,:)
    integer                     :: k
    integer                     :: k_chk
    type(evol_model_t), pointer :: em

    ! Read data from the MESA-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header and determine the version

    read(unit, *) n, global_data, version

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'File version', version/100._WP
120    format(3X,A,1X,F4.2,1X,A)
    endif

    ! Read the data

    select case (version)
    case (1)
       backspace(unit)
       n_cols = N_COLS_V0_01
    case (19)
       n_cols = N_COLS_V0_19
    case (100,101)
       n_cols = N_COLS_V1_0X
    case default
       $ABORT(Unrecognized MESA file version)
    end select

    allocate(point_data(n_cols,n))

    read_loop : do k = 1,n
       read(unit, *) k_chk, point_data(:,k)
       $ASSERT(k == k_chk,Index mismatch)
    end do read_loop

    close(unit)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', n, 'points'
130    format(3X,A,1X,I0,1X,A)
    endif

    ! Initialize the model

    call init_mesa_model(ml_p, global_data, point_data, version, em)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Return the pointer to the model

    ml => em

    ! Finish

    return

  end subroutine read_mesa_model

  !****

  subroutine init_mesa_model (ml_p, global_data, point_data, version, em)

    type(model_par_t), intent(in)            :: ml_p
    real(WP), intent(in)                     :: global_data(:)
    real(WP), intent(in)                     :: point_data(:,:)
    integer, intent(in)                      :: version
    type(evol_model_t), pointer, intent(out) :: em

    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    integer                     :: n
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: M_r(:)
    real(WP), allocatable       :: L_r(:)
    real(WP), allocatable       :: P(:)
    real(WP), allocatable       :: rho(:)
    real(WP), allocatable       :: T(:)
    real(WP), allocatable       :: N2(:)
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: eps(:)
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
    real(WP), allocatable       :: kap(:)
    real(WP), allocatable       :: kap_rho(:)
    real(WP), allocatable       :: kap_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: c_lum(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thn(:)
    real(WP), allocatable       :: c_thk(:)
    real(WP), allocatable       :: c_eps(:)
    real(WP), allocatable       :: int_rhoT(:)
    real(WP), allocatable       :: f_luan_t(:)
    real(WP), allocatable       :: f_luan_c(:)

    ! Extract data from the global and point arrays

    M_star = global_data(1)
    R_star = global_data(2)
    L_star = global_data(3)

    n = SIZE(point_data, 2)

    select case (version)
    case (1)
       call extract_data_v0_01_()
    case (19)
       call extract_data_v0_19_()
    case (100,101)
       call extract_data_v1_0X_()
    case default
       $ABORT(Unrecognized MESA memory version)
    end select

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, M_r)

    ! Calculate dimensionless structure data

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))
    allocate(c_lum(n))
    allocate(f_luan_t(n))
    
    where (x /= 0._WP)
       V_2 = G_GRAVITY*M_r*rho/(P*r*x**2)
       As = r**3*N2/(G_GRAVITY*M_r)
       U = 4._WP*PI*rho*r**3/M_r
       c_1 = (r/R_star)**3/(M_r/M_star)
       c_lum = (L_r/L_star)/x**3
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
       c_lum = 4._WP*PI*rho(1)*eps(1)*R_star**3/L_star
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kap*rho*L_star)
    c_thn = c_P*SQRT(G_GRAVITY*M_star/R_star**3)/(A_RADIATION*C_LIGHT*kap*T**3)
    c_thk = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star

    select case (version)
    case (101)
       c_eps = 4._WP*PI*rho*eps*R_star**3/L_star
    case default
       ! Note: technically incorrect because eps in versions < 1.01 includes eps_grav
       c_eps = 4._WP*PI*rho*eps*R_star**3/L_star
    end select

    int_rhoT = -integral(r(n:1:-1), rho(n:1:-1)*K_BOLTZMANN*T(n:1:-1)/M_PROTON)
    int_rhoT = int_rhoT(n:1:-1)
    where (x /= 0._WP)
       f_luan_t = 4._WP*PI*r**2*int_rhoT/L_r*SQRT(G_GRAVITY*M_star/R_star**3)
    elsewhere
       f_luan_t = 0._WP
    end where

    f_luan_c = c_P*M_PROTON/K_BOLTZMANN

    if (ml_p%uniform_rot) then
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

    call em%define(I_C_LUM, c_lum)
    call em%define(I_C_RAD, c_rad)
    call em%define(I_C_THN, c_thn)
    call em%define(I_C_THK, c_thk)
    call em%define(I_C_EPS, c_eps)

    call em%define(I_EPS_RHO, eps_rho)
    call em%define(I_EPS_T, eps_T)
       
    call em%define(I_KAP_RHO, kap_rho)
    call em%define(I_KAP_T, kap_T)

    call em%define(I_F_LUAN_T, f_luan_t)
    call em%define(I_F_LUAN_C, f_luan_c)

    call em%define(I_OMEGA_ROT, Omega_rot)

    ! Finish

    return

  contains

    subroutine extract_data_v0_01_ ()

      real(WP), allocatable :: eps_eps_rho(:)
      real(WP), allocatable :: eps_eps_T(:)
      integer               :: k

      $CHECK_BOUNDS(SIZE(point_data, 1),N_COLS_V0_01)

      ! Extract data from the version-0.01 point array

      r = point_data(1,:)
      M_r = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      L_r = point_data(3,:)
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(12,:)*point_data(10,:)/point_data(9,:)
      delta = point_data(11,:)/point_data(12,:)
      kap = point_data(13,:)
      kap_T = point_data(14,:)
      kap_rho = point_data(15,:)
      eps = point_data(16,:)
      eps_eps_T = point_data(17,:)
      eps_eps_rho = point_data(18,:)

      nabla_ad = p*delta/(rho*T*point_data(10,:))

      allocate(Omega_rot(n))
      Omega_rot = 0._WP

      ! Evaluate eps_rho and eps_T from eps_eps_*

      k = MAXLOC(ABS(eps_eps_T), DIM=1)

      if (ABS(eps_eps_T(k)) >= 1E-3_WP*ABS(eps(k))) then

         ! Note: technically incorrect because eps in versions < 1.01 includes eps_grav

         eps_rho = eps_eps_rho/eps
         eps_T = eps_eps_T/eps

         if(check_log_level('INFO')) then
            write(OUTPUT_UNIT, 120) 'Rescaled epsilon derivatives'
120         format(3X,A)
         endif

      else

         eps_rho = eps_eps_rho
         eps_T = eps_eps_T

      endif

      ! Finish

      return

    end subroutine extract_data_v0_01_

    subroutine extract_data_v0_19_ ()

      real(WP), allocatable :: eps_eps_rho(:)
      real(WP), allocatable :: eps_eps_T(:)

      $CHECK_BOUNDS(SIZE(point_data, 1),N_COLS_V0_19)

      ! Extract data from the version-0.19 point array

      r = point_data(1,:)
      M_r = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      L_r = point_data(3,:)
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kap = point_data(12,:)
      kap_T = point_data(13,:)
      kap_rho = point_data(14,:)
      eps = point_data(15,:)
      eps_eps_T = point_data(16,:)
      eps_eps_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Note: technically incorrect because eps in versions < 1.01 includes eps_grav

      allocate(eps_rho(n))
      allocate(eps_T(n))

      where (eps /= 0._WP)
         eps_rho = eps_eps_rho/eps
         eps_T = eps_eps_T/eps
      elsewhere
         eps_rho = 0._WP
         eps_T = 0._WP
      end where

      ! Finish

      return

    end subroutine extract_data_v0_19_

    subroutine extract_data_v1_0X_ ()

      real(WP), allocatable :: eps_eps_rho(:)
      real(WP), allocatable :: eps_eps_T(:)
      real(WP), allocatable :: kap_kap_T(:)
      real(WP), allocatable :: kap_kap_rho(:)

      $CHECK_BOUNDS(SIZE(point_data, 1),N_COLS_V1_0X)

      ! Extract data from the version-1.0X point array

      r = point_data(1,:)
      M_r = point_data(2,:)
      L_r = point_data(3,:)
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kap = point_data(12,:)
      kap_kap_T = point_data(13,:)
      kap_kap_rho = point_data(14,:)
      eps = point_data(15,:)
      eps_eps_T = point_data(16,:)
      eps_eps_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Note: technically incorrect in version 1.00, because eps in versions < 1.01 includes eps_grav

      allocate(eps_rho(n))
      allocate(eps_T(n))

      where (eps /= 0._WP)
         eps_rho = eps_eps_rho/eps
         eps_T = eps_eps_T/eps
      elsewhere
         eps_rho = 0._WP
         eps_T = 0._WP
      end where

      kap_T = kap_kap_T/kap
      kap_rho = kap_kap_rho/kap

      ! Finish

      return

    end subroutine extract_data_v1_0X_

  end subroutine init_mesa_model

end module gyre_mesa_file
