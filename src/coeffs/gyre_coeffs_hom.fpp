! Module   : gyre_coeffs_hom
! Purpose  : base structure coefficients for homogeneous compressible models
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

module gyre_coeffs_hom

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_coeffs
  use gyre_cocache

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1
    procedure :: ${NAME}_v
  $endsub

  type, extends(coeffs_t) :: coeffs_hom_t
     private
     real(WP) :: dt_Gamma_1
   contains
     private
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(tau_thm)
     $PROC_DECL(Omega_rot)
     procedure, public :: pi_c
     procedure, public :: is_zero
     procedure, public :: attach_cache
     procedure, public :: detach_cache
     procedure, public :: fill_cache
  end type coeffs_hom_t

  ! Interfaces

  interface coeffs_hom_t
     module procedure init_cf
  end interface coeffs_hom_t

  $if ($MPI)
  interface bcast
     module procedure bcast_cf
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: coeffs_hom_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function init_cf (Gamma_1) result (cf)

    real(WP), intent(in) :: Gamma_1
    type(coeffs_hom_t)   :: cf

    ! Construct the coeffs_hom

    cf%dt_Gamma_1 = Gamma_1

    ! Finish

    return

  end function init_cf

!****

  $if($MPI)

  subroutine bcast_cf (cf, root_rank)

    class(coeffs_hom_t), intent(inout) :: cf
    integer, intent(in)                :: root_rank

    ! Broadcast the coeffs_hom

    call bcast(cf%dt_Gamma_1, root_rank)

    ! Finish

    return

  end subroutine bcast_cf

  $endif

!****

  function V_1 (this, x) result (V)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: V

    ! Calculate V

    V = 2._WP*x**2/(1._WP - x**2)

    ! Finish

    return

  end function V_1

!****
  
  function V_v (this, x) result (V)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: V(SIZE(x))

    integer :: i

    ! Calculate V

    x_loop : do i = 1,SIZE(x)
       V(i) = this%V(x(i))
    end do x_loop

    ! Finish

    return

  end function V_v

!****

  function As_1 (this, x) result (As)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: As

    ! Calculate As

    As = -this%V(x)/this%dt_Gamma_1

    ! Finish

    return

  end function As_1

!****
  
  function As_v (this, x) result (As)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: As(SIZE(x))

    integer :: i

    ! Calculate As

    x_loop : do i = 1,SIZE(x)
       As(i) = this%As(x(i))
    end do x_loop

    ! Finish

    return

  end function As_v

!****

  function U_1 (this, x) result (U)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: U

    ! Calculate U

    U = 3._WP

    ! Finish

    return

  end function U_1

!****
  
  function U_v (this, x) result (U)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: U(SIZE(x))

    integer :: i

    ! Calculate U

    x_loop : do i = 1,SIZE(x)
       U(i) = this%U(x(i))
    end do x_loop

    ! Finish

    return

  end function U_v

!****

  function c_1_1 (this, x) result (c_1)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: c_1

    ! Calculate c_1

    c_1 = 1._WP

    ! Finish

    return

  end function c_1_1

!****
  
  function c_1_v (this, x) result (c_1)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: c_1(SIZE(x))

    integer :: i

    ! Calculate c_1

    x_loop : do i = 1,SIZE(x)
       c_1(i) = this%c_1(x(i))
    end do x_loop

    ! Finish

    return

  end function c_1_v

!****

  function Gamma_1_1 (this, x) result (Gamma_1)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: Gamma_1

    ! Calculate Gamma_1

    Gamma_1 = this%dt_Gamma_1

    ! Finish

    return

  end function Gamma_1_1

!****
  
  function Gamma_1_v (this, x) result (Gamma_1)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: Gamma_1(SIZE(x))

    integer :: i

    ! Calculate Gamma_1
    
    x_loop : do i = 1,SIZE(x)
       Gamma_1(i) = this%Gamma_1(x(i))
    end do x_loop

    ! Finish

    return

  end function Gamma_1_v

!****

  function nabla_ad_1 (this, x) result (nabla_ad)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: nabla_ad

    ! Calculate nabla_ad (assume ideal gas)

    nabla_ad = 2._WP/5._WP

    ! Finish

    return

  end function nabla_ad_1

!****
  
  function nabla_ad_v (this, x) result (nabla_ad)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: nabla_ad(SIZE(x))

    integer :: i

    ! Calculate nabla_ad
    
    x_loop : do i = 1,SIZE(x)
       nabla_ad(i) = this%nabla_ad(x(i))
    end do x_loop

    ! Finish

    return

  end function nabla_ad_v

!****

  function delta_1 (this, x) result (delta)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: delta

    ! Calculate delta (assume ideal gas)

    delta = 1._WP

    ! Finish

    return

  end function delta_1

!****
  
  function delta_v (this, x) result (delta)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: delta(SIZE(x))

    integer :: i

    ! Calculate delta
    
    x_loop : do i = 1,SIZE(x)
       delta(i) = this%delta(x(i))
    end do x_loop

    ! Finish

    return

  end function delta_v

!****

  $define $PROC $sub

  $local $NAME $1

  function ${NAME}_1 (this, x) result ($NAME)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: $NAME

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_1

!****

  function ${NAME}_v (this, x) result ($NAME)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: $NAME(SIZE(x))

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_v

  $endsub

  $PROC(c_rad)
  $PROC(dc_rad)
  $PROC(c_thm)
  $PROC(c_dif)
  $PROC(c_eps_ad)
  $PROC(c_eps_S)
  $PROC(nabla)
  $PROC(kappa_S)
  $PROC(kappa_ad)
  $PROC(tau_thm)

!****

  function Omega_rot_1 (this, x) result (Omega_rot)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: Omega_rot

    ! Calculate Omega_rot (no rotation)

    Omega_rot = 0._WP

    ! Finish

    return

  end function Omega_rot_1

!****
  
  function Omega_rot_v (this, x) result (Omega_rot)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: Omega_rot(SIZE(x))

    integer :: i

    ! Calculate Omega_rot
    
    x_loop : do i = 1,SIZE(x)
       Omega_rot(i) = this%Omega_rot(x(i))
    end do x_loop

    ! Finish

    return

  end function Omega_rot_v

!****

  function pi_c (this)

    class(coeffs_hom_t), intent(in) :: this
    real(WP)                        :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c = 2._WP

    ! Finish

    return

  end function pi_c

!****

  function is_zero (this, x)

    class(coeffs_hom_t), intent(in) :: this
    real(WP), intent(in)            :: x
    logical                         :: is_zero

    ! Determine whether the point at x has a vanishing pressure and/or density

    is_zero = x >= 1._WP

    ! Finish

    return

  end function is_zero

!****

  subroutine attach_cache (this, cc)

    class(coeffs_hom_t), intent(inout)    :: this
    class(cocache_t), pointer, intent(in) :: cc

    ! Attach a coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine attach_cache

!****

  subroutine detach_cache (this)

    class(coeffs_hom_t), intent(inout) :: this

    ! Detach the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine detach_cache

!****

  subroutine fill_cache (this, x)

    class(coeffs_hom_t), intent(inout) :: this
    real(WP), intent(in)               :: x(:)

    ! Fill the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine fill_cache

end module gyre_coeffs_hom
