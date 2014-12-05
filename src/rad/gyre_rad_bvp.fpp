! Module   : gyre_rad_bvp
! Purpose  : boundary-value solver (adiabatic radial)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_rad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_model
  use gyre_mode
  use gyre_ext
  use gyre_ivp
  use gyre_sysmtx
  use gyre_modepar
  use gyre_oscpar
  use gyre_numpar
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(r_bvp_t) :: rad_bvp_t
   contains
     private
     procedure, public :: recon => recon_
  end type rad_bvp_t

  ! Interfaces

  interface rad_bvp_t
     module procedure rad_bvp_t_
  end interface rad_bvp_t

  ! Access specifiers

  private

  public :: rad_bvp_t

  ! Procedures

contains

  function rad_bvp_t_ (x, ml, mp, op, np) result (bp)

    use gyre_rad_jacob
    use gyre_rad_bound

    use gyre_magnus_ivp
    use gyre_colloc_ivp
    use gyre_findiff_ivp

    use gyre_block_sysmtx
 
    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    type(numpar_t), intent(in)          :: np
    type(rad_bvp_t), target             :: bp

    type(rad_jacob_t)              :: jc
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(rad_bound_t)              :: bd
    class(r_ivp_t), allocatable    :: iv
    class(r_sysmtx_t), allocatable :: sm
 
    ! Construct the rad_bvp_t

    ! Initialize the jacobian

    select case (op%variables_type)
    case ('DZIEM')
       jc = rad_jacob_t(ml, mp, 'DZIEM')
    case ('JCD')
       jc = rad_jacob_t(ml, mp, 'JCD')
    case ('MIX')
       jc = rad_jacob_t(ml, mp, 'MIX')
    case default
       $ABORT(Invalid variables_type)
    end select

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    print *,'Set:',x_i,x_o

    select case (op%outer_bound_type)
    case ('ZERO')
       bd = rad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'ZERO')
    case ('DZIEM')
       bd = rad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'DZIEM')
    case ('UNNO')
       bd = rad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'UNNO')
    case ('JCD')
       bd = rad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'JCD')
    case default
       $ABORT(Invalid bound_type)
    end select

    ! Initialize the IVP solver

    select case (np%ivp_solver_type)
    case ('MAGNUS_GL2')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL6'))
    case ('COLLOC_GL2')
       allocate(iv, SOURCE=r_colloc_ivp_t(jc, 'GL2'))
    case ('COLLOC_GL4')
       allocate(iv, SOURCE=r_colloc_ivp_t(jc, 'GL4'))
    case ('FINDIFF')
       allocate(iv, SOURCE=r_findiff_ivp_t(jc))
    case default
       $ABORT(Invalid ivp_solver_type)
    end select

    ! Initialize the system matrix

    if (np%use_banded) then
       $ABORT(Not yet implemented)
    else
       allocate(sm, SOURCE=r_block_sysmtx_t(n-1, jc%n_e, bd%n_i, bd%n_o))
    endif

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t(x, ml, jc, bd, iv, sm)

    ! Finish

    return

  end function rad_bvp_t_

!****

  subroutine recon_ (this, omega, x, x_ref, y, y_ref, discrim)

    class(rad_bvp_t), intent(inout) :: this
    real(WP), intent(in)            :: omega
    real(WP), intent(in)            :: x(:)
    real(WP), intent(in)            :: x_ref
    real(WP), intent(out)           :: y(:,:)
    real(WP), intent(out)           :: y_ref(:)
    type(r_ext_t), intent(out)      :: discrim

    real(WP) :: y_(2,SIZE(x))
    real(WP) :: y_ref_(2)
    integer  :: n
    integer  :: i
    real(WP) :: y_4_x(SIZE(x))
    real(WP) :: eul_phi(SIZE(x))

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    $CHECK_BOUNDS(SIZE(y_ref),6)

    ! Reconstruct the solution

    call this%r_bvp_t%recon(omega, x, x_ref, y_, y_ref_, discrim)

    ! Convert to the canonical (6-variable) solution

    n = SIZE(x)

    !$OMP PARALLEL DO 
    do i = 1, n
       y(1:2,i) = MATMUL(this%jc%T(x(i), omega, .TRUE.), y_(:,i))
       y(4,i) = -y(1,i)*this%ml%U(x(i))
       y(5:6,i) = 0._WP
    end do

    ! Reconstruct the potential by integrating the gravity

    where (x /= 0._WP)
       y_4_x = y(4,:)/x
    elsewhere
       y_4_x = 0._WP
    end where

    eul_phi = integral(x, y_4_x/this%ml%c_1(x))

    y(3,:) = this%ml%c_1(x)*(eul_phi - eul_phi(n))
    y(2,:) = y(2,:) + y(3,:)

    ! Note: y_ref(3) is set to zero; need to fix this

    y_ref(1:2) = MATMUL(this%jc%T(x_ref, omega, .TRUE.), y_ref_)
    y_ref(3) = 0._WP
    y_ref(4) = -y_ref_(1)*this%ml%U(x_ref)
    y_ref(5:6) = 0._WP

    ! Finish

    return

  end subroutine recon_

! ! !****

! !   function mode_c_ (this, omega) result (md)

! !     class(rad_bvp_t), target, intent(inout) :: this
! !     complex(WP), intent(in)                 :: omega
! !     type(mode_t)                            :: md

! !     real(WP), allocatable    :: x(:)
! !     complex(WP), allocatable :: y(:,:)
! !     real(WP)                 :: x_ref
! !     complex(WP)              :: y_ref(this%n_e)
! !     type(ext_complex_t)      :: discrim
! !     integer                  :: n
! !     integer                  :: i
! !     complex(WP), allocatable :: y_c(:,:)
! !     complex(WP), allocatable :: y_4_x(:)
! !     complex(WP), allocatable :: eul_phi(:)
! !     complex(WP)              :: y_c_ref(6)

! !     ! Reconstruct the solution

! !     call this%recon_(omega, x, y, x_ref, y_ref, discrim, .FALSE.)

! !     ! Calculate canonical variables

! !     n = SIZE(x)

! !     allocate(y_c(6,n))

! !     !$OMP PARALLEL DO 
! !     do i = 1,n
! !        y_c(1:2,i) = MATMUL(this%jc%trans_matrix(x(i), omega, .TRUE.), y(:,i))
! !        y_c(4,i) = -y_c(1,i)*this%ml%U(x(i))
! !        y_c(5:6,i) = 0._WP
! !     end do

! !     allocate(y_4_x(n))

! !     where (x /= 0._WP)
! !        y_4_x = y_c(4,:)/x
! !     elsewhere
! !        y_4_x = 0._WP
! !     end where

! !     eul_phi = integral(x, y_4_x/this%ml%c_1(x))

! !     y_c(3,:) = this%ml%c_1(x)*(eul_phi - eul_phi(n))
! !     y_c(2,:) = y_c(2,:) + y_c(3,:)

! !     y_c_ref(1:2) = MATMUL(this%jc%trans_matrix(x_ref, omega, .TRUE.), y_ref)
! !     y_c_ref(3) = 0._WP
! !     y_c_ref(4) = -y_c_ref(1)*this%ml%U(x_ref)
! !     y_c_ref(5:6) = 0._WP

! !     ! Initialize the mode

! !     md = mode_t(this%ml, this%mp, this%op, omega, discrim, x, y_c, x_ref, y_c_ref)

! !     ! Finish

! !     return

! !   end function mode_c_

! !****

!   subroutine build_ (this, omega)

!     class(rad_bvp_t), target, intent(inout) :: this
!     complex(WP), intent(in)                 :: omega

!     ! Set up the sysmtx

!     call this%ml%attach_cache(this%cc)

!     call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega), ext_complex_t(1._WP))
!     call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega), ext_complex_t(1._WP))

!     call this%sh%shoot(omega, this%x, this%sm)

!     call this%ml%detach_cache()

!     call this%sm%scale_rows()

!     ! Finish

!     return

!   end subroutine build_

! !****

!   subroutine recon_ (this, omega, x, y, x_ref, y_ref, discrim, use_real)

!     class(rad_bvp_t), intent(inout)       :: this
!     complex(WP), intent(in)               :: omega
!     real(WP), allocatable, intent(out)    :: x(:)
!     complex(WP), allocatable, intent(out) :: y(:,:)
!     real(WP), intent(out)                 :: x_ref
!     complex(WP), intent(out)              :: y_ref(:)
!     type(ext_complex_t), intent(out)      :: discrim
!     logical, optional, intent(in)         :: use_real

!     complex(WP) :: b(this%n_e*this%n)
!     complex(WP) :: y_sh(this%n_e,this%n)
!     logical     :: same_grid
!     complex(WP) :: y_ref_(this%n_e,1)

!     $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

!     ! Reconstruct the solution on the shooting grid

!     call this%build_(omega)

!     call this%sm%null_vector(b, discrim, use_real, this%np%use_banded)

!     y_sh = RESHAPE(b, SHAPE(y_sh))

!     ! Build the recon grid

!     this%recon_gp%omega_a = REAL(omega)
!     this%recon_gp%omega_b = REAL(omega)

!     call build_grid(this%recon_gp, this%ml, this%mp, this%x, x)

!     if(SIZE(x) == SIZE(this%x)) then
!        same_grid = ALL(x == this%x)
!     else
!        same_grid = .FALSE.
!     endif

!     ! Reconstruct the full solution

!     if(same_grid) then

!        y = y_sh

!     else

!        allocate(y(this%n_e,SIZE(x)))

!        call this%sh%recon(omega, this%x, y_sh, x, y)

!     endif

!     ! Reconstruct the solution at x_ref
    
!     x_ref = MIN(MAX(this%op%x_ref, this%x(1)), this%x(this%n))

!     call this%sh%recon(omega, this%x, y_sh, [x_ref], y_ref_)

!     y_ref = y_ref_(:,1)

!     ! Finish

!     return

!   end subroutine recon_

! !****

!   function model_ (this) result (ml)

!     class(rad_bvp_t), intent(in) :: this
!     class(model_t), pointer      :: ml

!     ! Return the model pointer

!     ml => this%ml

!     ! Finish

!     return

!   end function model_

! !****

!   function mode_r_ (this, omega) result (md)

!     class(rad_bvp_t), target, intent(inout) :: this
!     real(WP), intent(in)                    :: omega
!     type(mode_t)                            :: md

!     complex(WP)              :: omega_c
!     real(WP), allocatable    :: x(:)
!     complex(WP), allocatable :: y(:,:)
!     real(WP)                 :: x_ref
!     complex(WP)              :: y_ref(this%n_e)
!     type(ext_complex_t)      :: discrim_c
!     integer                  :: n
!     integer                  :: i
!     complex(WP), allocatable :: y_c(:,:)
!     complex(WP), allocatable :: y_4_x(:)
!     complex(WP), allocatable :: eul_phi(:)
!     complex(WP)              :: y_c_ref(6)

!     ! Reconstruct the solution

!     omega_c = CMPLX(omega, KIND=WP)

!     call this%recon_(omega_c, x, y, x_ref, y_ref, discrim_c, .TRUE.)

!     ! Calculate canonical variables

!     n = SIZE(x)

!     allocate(y_c(6,n))

!     !$OMP PARALLEL DO 
!     do i = 1,n
!        y_c(1:2,i) = MATMUL(this%jc%trans_matrix(x(i), omega_c, .TRUE.), y(:,i))
!        y_c(4,i) = -y_c(1,i)*this%ml%U(x(i))
!        y_c(5:6,i) = 0._WP
!     end do

!     allocate(y_4_x(n))

!     where (x /= 0._WP)
!        y_4_x = y_c(4,:)/x
!     elsewhere
!        y_4_x = 0._WP
!     end where

!     eul_phi = integral(x, y_4_x/this%ml%c_1(x))

!     y_c(3,:) = this%ml%c_1(x)*(eul_phi - eul_phi(n))
!     y_c(2,:) = y_c(2,:) + y_c(3,:)

!     y_c_ref(1:2) = MATMUL(this%jc%trans_matrix(x_ref, omega_c, .TRUE.), y_ref)
!     y_c_ref(3) = 0._WP
!     y_c_ref(4) = -y_c_ref(1)*this%ml%U(x_ref)
!     y_c_ref(5:6) = 0._WP

!     ! Initialize the mode

!     md = mode_t(this%ml, this%mp, this%op, omega_c, discrim_c, x, y_c, x_ref, y_c_ref)

!     ! Finish

!     return

!   end function mode_r_

! !****

!   function mode_c_ (this, omega) result (md)

!     class(rad_bvp_t), target, intent(inout) :: this
!     complex(WP), intent(in)                 :: omega
!     type(mode_t)                            :: md

!     real(WP), allocatable    :: x(:)
!     complex(WP), allocatable :: y(:,:)
!     real(WP)                 :: x_ref
!     complex(WP)              :: y_ref(this%n_e)
!     type(ext_complex_t)      :: discrim
!     integer                  :: n
!     integer                  :: i
!     complex(WP), allocatable :: y_c(:,:)
!     complex(WP), allocatable :: y_4_x(:)
!     complex(WP), allocatable :: eul_phi(:)
!     complex(WP)              :: y_c_ref(6)

!     ! Reconstruct the solution

!     call this%recon_(omega, x, y, x_ref, y_ref, discrim, .FALSE.)

!     ! Calculate canonical variables

!     n = SIZE(x)

!     allocate(y_c(6,n))

!     !$OMP PARALLEL DO 
!     do i = 1,n
!        y_c(1:2,i) = MATMUL(this%jc%trans_matrix(x(i), omega, .TRUE.), y(:,i))
!        y_c(4,i) = -y_c(1,i)*this%ml%U(x(i))
!        y_c(5:6,i) = 0._WP
!     end do

!     allocate(y_4_x(n))

!     where (x /= 0._WP)
!        y_4_x = y_c(4,:)/x
!     elsewhere
!        y_4_x = 0._WP
!     end where

!     eul_phi = integral(x, y_4_x/this%ml%c_1(x))

!     y_c(3,:) = this%ml%c_1(x)*(eul_phi - eul_phi(n))
!     y_c(2,:) = y_c(2,:) + y_c(3,:)

!     y_c_ref(1:2) = MATMUL(this%jc%trans_matrix(x_ref, omega, .TRUE.), y_ref)
!     y_c_ref(3) = 0._WP
!     y_c_ref(4) = -y_c_ref(1)*this%ml%U(x_ref)
!     y_c_ref(5:6) = 0._WP

!     ! Initialize the mode

!     md = mode_t(this%ml, this%mp, this%op, omega, discrim, x, y_c, x_ref, y_c_ref)

!     ! Finish

!     return

!   end function mode_c_

! !****

!   subroutine build_ (this, omega)

!     class(rad_bvp_t), target, intent(inout) :: this
!     complex(WP), intent(in)                 :: omega

!     ! Set up the sysmtx

!     call this%ml%attach_cache(this%cc)

!     call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega), ext_complex_t(1._WP))
!     call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega), ext_complex_t(1._WP))

!     call this%sh%shoot(omega, this%x, this%sm)

!     call this%ml%detach_cache()

!     call this%sm%scale_rows()

!     ! Finish

!     return

!   end subroutine build_

! !****

!   subroutine recon_ (this, omega, x, y, x_ref, y_ref, discrim, use_real)

!     class(rad_bvp_t), intent(inout)       :: this
!     complex(WP), intent(in)               :: omega
!     real(WP), allocatable, intent(out)    :: x(:)
!     complex(WP), allocatable, intent(out) :: y(:,:)
!     real(WP), intent(out)                 :: x_ref
!     complex(WP), intent(out)              :: y_ref(:)
!     type(ext_complex_t), intent(out)      :: discrim
!     logical, optional, intent(in)         :: use_real

!     complex(WP) :: b(this%n_e*this%n)
!     complex(WP) :: y_sh(this%n_e,this%n)
!     logical     :: same_grid
!     complex(WP) :: y_ref_(this%n_e,1)

!     $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

!     ! Reconstruct the solution on the shooting grid

!     call this%build_(omega)

!     call this%sm%null_vector(b, discrim, use_real, this%np%use_banded)

!     y_sh = RESHAPE(b, SHAPE(y_sh))

!     ! Build the recon grid

!     this%recon_gp%omega_a = REAL(omega)
!     this%recon_gp%omega_b = REAL(omega)

!     call build_grid(this%recon_gp, this%ml, this%mp, this%x, x)

!     if(SIZE(x) == SIZE(this%x)) then
!        same_grid = ALL(x == this%x)
!     else
!        same_grid = .FALSE.
!     endif

!     ! Reconstruct the full solution

!     if(same_grid) then

!        y = y_sh

!     else

!        allocate(y(this%n_e,SIZE(x)))

!        call this%sh%recon(omega, this%x, y_sh, x, y)

!     endif

!     ! Reconstruct the solution at x_ref
    
!     x_ref = MIN(MAX(this%op%x_ref, this%x(1)), this%x(this%n))

!     call this%sh%recon(omega, this%x, y_sh, [x_ref], y_ref_)

!     y_ref = y_ref_(:,1)

!     ! Finish

!     return

!   end subroutine recon_

! !****

!   function model_ (this) result (ml)

!     class(rad_bvp_t), intent(in) :: this
!     class(model_t), pointer      :: ml

!     ! Return the model pointer

!     ml => this%ml

!     ! Finish

!     return

!   end function model_

end module gyre_rad_bvp
