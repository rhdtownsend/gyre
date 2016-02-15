! Module   : gyre_mesa_file
! Purpose  : read MESA files
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

module gyre_mesa_file

  ! Uses

  use core_kinds
  use core_order

  use gyre_constants
  use gyre_evol_model
  use gyre_model
  use gyre_model_par
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_mesa_model
  public :: read_mesa_data

  ! Procedures

contains

  subroutine read_mesa_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: P(:)
    real(WP), allocatable       :: rho(:)
    real(WP), allocatable       :: T(:)
    real(WP), allocatable       :: N2(:)
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kappa(:)
    real(WP), allocatable       :: kappa_rho(:)
    real(WP), allocatable       :: kappa_T(:)
    real(WP), allocatable       :: epsilon(:)
    real(WP), allocatable       :: epsilon_rho(:)
    real(WP), allocatable       :: epsilon_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    integer                     :: n
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thm(:)
    real(WP), allocatable       :: c_dif(:)
    real(WP), allocatable       :: c_eps_ad(:)
    real(WP), allocatable       :: c_eps_S(:)
    real(WP), allocatable       :: kappa_ad(:)
    real(WP), allocatable       :: kappa_S(:)
    type(evol_model_t), pointer :: em

    ! Read data from the MESA-format file

    call read_mesa_data(ml_p%file, M_star, R_star, L_star, r, m, P, rho, T, &
                        N2, Gamma_1, nabla_ad, delta, nabla,  &
                        kappa, kappa_rho, kappa_T, &
                        epsilon, epsilon_rho, epsilon_T, &
                        Omega_rot)

    ! Calculate dimensionless structure data

    x = r/R_star

    n = SIZE(x)

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(P*r*x**2)
       As = r**3*N2/(G_GRAVITY*m)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

    kappa_ad = nabla_ad*kappa_T + kappa_rho/Gamma_1
    kappa_S = kappa_T - delta*kappa_rho

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kappa*rho*L_star)
    c_thm = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_dif = (kappa_ad-4._WP*nabla_ad)*V_2*x**2*nabla + V_2*x**2*nabla_ad

    c_eps_ad = 4._WP*PI*rho*(nabla_ad*epsilon_T + epsilon_rho/Gamma_1)*R_star**3/L_star
    c_eps_S = 4._WP*PI*rho*(epsilon_T - delta*epsilon_rho)*R_star**3/L_star

    if (ml_p%uniform_rot) then
       Omega_rot = ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%set_V_2(V_2)
    call em%set_As(As)
    call em%set_U(U)
    call em%set_c_1(c_1)

    call em%set_Gamma_1(Gamma_1)
    call em%set_delta(delta)
    call em%set_nabla_ad(nabla_ad)
    call em%set_nabla(nabla)
    call em%set_beta_rad(beta_rad)

    call em%set_c_rad(c_rad)
    call em%set_c_thm(c_thm)
    call em%set_c_dif(c_dif)
    call em%set_c_eps_ad(c_eps_ad)
    call em%set_c_eps_S(c_eps_S)
    call em%set_kappa_ad(kappa_ad)
    call em%set_kappa_S(kappa_S)

    call em%set_Omega_rot(Omega_rot)

    ! Return a pointer

    ml => em

    ! Finish

    return

  end subroutine read_mesa_model

  !****

  subroutine read_mesa_data (file, M_star, R_star, L_star, r, m, p, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, &
                             Omega_rot)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: M_star
    real(WP), intent(out)              :: R_star
    real(WP), intent(out)              :: L_star
    real(WP), allocatable, intent(out) :: r(:)
    real(WP), allocatable, intent(out) :: m(:)
    real(WP), allocatable, intent(out) :: p(:)
    real(WP), allocatable, intent(out) :: rho(:)
    real(WP), allocatable, intent(out) :: T(:)
    real(WP), allocatable, intent(out) :: N2(:)
    real(WP), allocatable, intent(out) :: Gamma_1(:)
    real(WP), allocatable, intent(out) :: nabla_ad(:)
    real(WP), allocatable, intent(out) :: delta(:)
    real(WP), allocatable, intent(out) :: nabla(:)
    real(WP), allocatable, intent(out) :: kappa(:)
    real(WP), allocatable, intent(out) :: kappa_rho(:)
    real(WP), allocatable, intent(out) :: kappa_T(:)
    real(WP), allocatable, intent(out) :: epsilon(:)
    real(WP), allocatable, intent(out) :: epsilon_rho(:)
    real(WP), allocatable, intent(out) :: epsilon_T(:)
    real(WP), allocatable, intent(out) :: Omega_rot(:)

    integer  :: unit
    integer  :: n
    integer  :: version
    
    ! Read data from the MESA-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file', TRIM(file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header to determine the version

    read(unit, *) n, M_star, R_star, L_star, version

    select case (version)
    case (1)

       backspace(unit)

       if (check_log_level('INFO')) then
         write(OUTPUT_UNIT, 100) 'Detected version 0.01 file'
      endif

       call read_mesa_data_v0_01_()

    case (19)

       if (check_log_level('INFO')) then
         write(OUTPUT_UNIT, 100) 'Detected version 0.19 file'
      endif

       call read_mesa_data_v0_19_()

    case (100)

      if (check_log_level('INFO')) then
         write(OUTPUT_UNIT, 100) 'Detected version 1.00 file'
      endif

       call read_mesa_data_v1_00_()

    case default

       $ABORT(Unrecognized MESA file version)

    end select

    close(unit)

    ! Finish

    return

  contains

    subroutine read_mesa_data_v0_01_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-0.01 file (this is the old variant
      ! with no version number; the 1 comes from reading the first
      ! field of the following record)

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      m = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(12,:)*point_data(10,:)/point_data(9,:)
      delta = point_data(11,:)/point_data(12,:)
      kappa = point_data(13,:)
      kappa_T = point_data(14,:)
      kappa_rho = point_data(15,:)
      epsilon = point_data(16,:)
      epsilon_T = point_data(17,:)
      epsilon_rho = point_data(18,:)

      nabla_ad = p*delta/(rho*T*point_data(10,:))

      allocate(Omega_rot(n))
      Omega_rot = 0._WP

      ! Decide whether epsilon_T and epsilon_rho need rescaling

      k = MAXLOC(ABS(epsilon_T), DIM=1)

      if (ABS(epsilon_T(k)) < 1E-3*ABS(epsilon(k))) then

         epsilon_T = epsilon_T*epsilon
         epsilon_rho = epsilon_rho*epsilon

         if(check_log_level('INFO')) then
            write(OUTPUT_UNIT, 120) 'Rescaled epsilon derivatives'
120         format(3X,A)
         endif

      endif

      ! Finish

      return

    end subroutine read_mesa_data_v0_01_

    subroutine read_mesa_data_v0_19_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-0.01 file (this is the old variant
      ! with no version number; the 19 comes from reading the column count)

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      m = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kappa = point_data(12,:)
      kappa_T = point_data(13,:)
      kappa_rho = point_data(14,:)
      epsilon = point_data(15,:)
      epsilon_T = point_data(16,:)
      epsilon_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_v0_19_

    subroutine read_mesa_data_v1_00_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-1.00 file

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      m = point_data(2,:)
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kappa = point_data(12,:)
      kappa_T = point_data(13,:)/point_data(12,:)
      kappa_rho = point_data(14,:)/point_data(12,:)
      epsilon = point_data(15,:)
      epsilon_T = point_data(16,:)
      epsilon_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_v1_00_

  end subroutine read_mesa_data

end module gyre_mesa_file
