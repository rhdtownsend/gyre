! Module   : gyre_mesa_file
! Purpose  : read MESA files
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

module gyre_mesa_file

  ! Uses

  use core_kinds
  use core_order

  use gyre_constants
  use gyre_model
  use gyre_evol_model
  use gyre_scons_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface read_mesa_model
     module procedure read_mesa_model_evol_
     module procedure read_mesa_model_scons_
  end interface read_mesa_model

  ! Access specifiers

  private

  public :: read_mesa_model
  public :: read_mesa_data

  ! Procedures

contains

  subroutine read_mesa_model_evol_ (file, deriv_type, regularize, ml, x)

    character(*), intent(in)                     :: file
    character(*), intent(in)                     :: deriv_type 
    logical, intent(in)                          :: regularize
    type(evol_model_t), intent(out)              :: ml
    real(WP), allocatable, optional, intent(out) :: x(:)

    integer               :: unit
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    integer               :: n_cols
    real(WP), allocatable :: var(:,:)
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: T(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: nabla_ad(:)
    real(WP), allocatable :: delta(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: Omega_rot(:)
    logical               :: add_center

    ! Read data from the MESA-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file', TRIM(file)
100    format(A,1X,A)
    endif

    call read_mesa_data(file, M_star, R_star, L_star, r, m, p, rho, T, &
                        N2, Gamma_1, nabla_ad, delta, nabla,  &
                        kappa, kappa_rho, kappa_T, &
                        epsilon, epsilon_rho, epsilon_T, &
                        Omega_rot)

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    if (check_log_level('INFO')) then
       if (regularize) write(OUTPUT_UNIT, 110) 'Regularizing'
       if (add_center) write(OUTPUT_UNIT, 110) 'Adding central point'
110    format(3X,A)
    endif

    ! Initialize the model

    ml = evol_model_t(M_star, R_star, L_star, r, m, p, rho, T, &
                      N2, Gamma_1, nabla_ad, delta, Omega_rot, &
                      nabla, kappa, kappa_rho, kappa_T, &
                      epsilon, epsilon_rho, epsilon_T, &
                      deriv_type, regularize=regularize, add_center=add_center)

    ! Set up the grid

    if(PRESENT(x)) then
       if(add_center) then
          x = [0._WP,r/R_star]
       else
          x = r/R_star
       endif
    endif
    ! Finish

    return

  end subroutine read_mesa_model_evol_

!****

  subroutine read_mesa_model_scons_ (file, ml, x)

    character(*), intent(in)                     :: file
    type(scons_model_t), intent(out)             :: ml
    real(WP), allocatable, optional, intent(out) :: x(:)

    integer               :: unit
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    integer               :: n_cols
    real(WP), allocatable :: var(:,:)
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: T(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: nabla_ad(:)
    real(WP), allocatable :: delta(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: Omega_rot(:)
    logical               :: add_center

    ! Read data from the MESA-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file', TRIM(file)
100    format(A,1X,A)
    endif

    call read_mesa_data(file, M_star, R_star, L_star, r, m, p, rho, T, &
                        N2, Gamma_1, nabla_ad, delta, nabla,  &
                        kappa, kappa_rho, kappa_T, &
                        epsilon, epsilon_rho, epsilon_T, &
                        Omega_rot)

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    if (check_log_level('INFO')) then
       if (add_center) write(OUTPUT_UNIT, 110) 'Adding central point'
110    format(3X,A)
    endif

    ! Initialize the model

    ml = scons_model_t(M_star, R_star, L_star, r, m, p, rho, T, &
                       N2, Gamma_1, nabla_ad, delta, Omega_rot, &
                       add_center=add_center)

    ! Set up the grid

    if(PRESENT(x)) then
       if(add_center) then
          x = [0._WP,r/R_star]
       else
          x = r/R_star
       endif
    endif
    ! Finish

    return

  end subroutine read_mesa_model_scons_

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
    integer  :: n_cols
    
    ! Read data from the MESA-format file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *) n, M_star, R_star, L_star, n_cols

    if(n_cols == 1) then

       ! Old variant (n_cols not specified)

       n_cols = 18

       backspace(unit)

       if(check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Detected old-variant file'
100       format(3X,A)
       endif

       call read_mesa_data_old_()

    else

       ! New variant (n_cols specified)

       if(check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Detected new-variant file'
       endif

       call read_mesa_data_new_()

    endif

    close(unit)

    ! Finish

    return

  contains

    subroutine read_mesa_data_old_ ()

      real(WP), allocatable :: var(:,:)
      integer               :: k
      integer               :: k_chk
      integer, allocatable  :: ind(:)

      ! Read data from the old-variant file

      allocate(var(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, var(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      close(unit)

      ind = unique_indices(var(1,:))

      if (SIZE(ind) < n) then

         if(check_log_level('WARN')) then
            write(OUTPUT_UNIT, 110) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
110         format('!!',1X,A)
         endif

         n = SIZE(var, 2)

      endif
       
      var = var(:,ind)

      r = var(1,:)
      m = var(2,:)/(1._WP+var(2,:))*M_star
      p = var(4,:)
      T = var(5,:)
      rho = var(6,:)
      nabla = var(7,:)
      N2 = var(8,:)
      Gamma_1 = var(12,:)*var(10,:)/var(9,:)
      delta = var(11,:)/var(12,:)
      kappa = var(13,:)
      kappa_T = var(14,:)
      kappa_rho = var(15,:)
      epsilon = var(16,:)
      epsilon_T = var(17,:)
      epsilon_rho = var(18,:)

      nabla_ad = p*delta/(rho*T*var(10,:))

      allocate(Omega_rot(n))
      Omega_rot = 0._WP

      ! Decide whether epsilon_T and epsilon_rho need rescaling

      k = MAXLOC(ABS(epsilon_T), DIM=1)

      if(ABS(epsilon_T(k)) < 1E-3*ABS(epsilon(k))) then

         epsilon_T = epsilon_T*epsilon
         epsilon_rho = epsilon_rho*epsilon

         if(check_log_level('INFO')) then
            write(OUTPUT_UNIT, 120) 'Rescaled epsilon derivatives'
120         format(3X,A)
         endif

      endif

      ! Finish

      return

    end subroutine read_mesa_data_old_

    subroutine read_mesa_data_new_ ()

      real(WP), allocatable :: var(:,:)
      integer               :: k
      integer               :: k_chk
      integer, allocatable  :: ind(:)

      $ASSERT(n_cols >= 18,Too few columns)

      ! Read data from the new-variant file

      allocate(var(n_cols-1,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, var(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      close(unit)

      ind = unique_indices(var(1,:))

      if (SIZE(ind) < n) then

         if(check_log_level('WARN')) then
            write(OUTPUT_UNIT, 110) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
110         format('!!',1X,A)
         endif

         n = SIZE(var, 2)

      endif
       
      var = var(:,ind)

      r = var(1,:)
      m = var(2,:)/(1._WP+var(2,:))*M_star
      p = var(4,:)
      T = var(5,:)
      rho = var(6,:)
      nabla = var(7,:)
      N2 = var(8,:)
      Gamma_1 = var(9,:)
      nabla_ad = var(10,:)
      delta = var(11,:)
      kappa = var(12,:)
      kappa_T = var(13,:)
      kappa_rho = var(14,:)
      epsilon = var(15,:)
      epsilon_T = var(16,:)
      epsilon_rho = var(17,:)
      Omega_rot = var(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_new_

  end subroutine read_mesa_data

end module gyre_mesa_file
