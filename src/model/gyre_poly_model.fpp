! Module   : gyre_poly_model
! Purpose  : stellar polytropic model
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

module gyre_poly_model

  ! Uses

  use core_kinds
  use core_parallel
  use core_spline

  use gyre_model
  use gyre_cocache

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  type, extends(model_t) :: poly_model_t
     private
     type(spline_t)   :: sp_Theta
     type(spline_t)   :: sp_dTheta
     type(spline_t)   :: sp_Omega_rot
     real(WP)         :: dt_Gamma_1
     real(WP), public :: n_poly
     real(WP), public :: xi_1
     real(WP), public :: Omega_uni
     logical, public  :: uniform_rot
   contains
     private
     $PROC_DECL(V_2)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(D)
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
     procedure, public :: pi_c => pi_c_
     procedure, public :: is_zero => is_zero_
     procedure, public :: attach_cache => attach_cache_
     procedure, public :: detach_cache => detach_cache_
     procedure, public :: fill_cache => fill_cache_
  end type poly_model_t

  ! Interfaces

  interface poly_model_t
     module procedure poly_model_t_
  end interface poly_model_t

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: poly_model_t
  $if ($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function poly_model_t_ (xi, Theta, dTheta, Omega_rot, n_poly, Gamma_1, deriv_type) result (ml)

    real(WP), intent(in)         :: xi(:)
    real(WP), intent(in)         :: Theta(:)
    real(WP), intent(in)         :: dTheta(:)
    real(WP), intent(in)         :: Omega_rot(:)
    real(WP), intent(in)         :: n_poly
    real(WP), intent(in)         :: Gamma_1
    character(LEN=*), intent(in) :: deriv_type
    type(poly_model_t)           :: ml

    integer  :: n
    real(WP) :: d2Theta(SIZE(xi))

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(Omega_rot),SIZE(xi))

    ! Construct the poly_model_t from the polytrope functions

    n = SIZE(xi)

    if(n_poly /= 0._WP) then

       where (xi /= 0._WP)
          d2Theta = -2._WP*dTheta/xi - Theta**n_poly
       elsewhere
          d2Theta = -1._WP/3._WP
       end where

    else

       d2Theta = -1._WP/3._WP

    endif

    ml%sp_Theta = spline_t(xi, Theta, dTheta)
    ml%sp_dTheta = spline_t(xi, dTheta, d2Theta)

    ml%sp_Omega_rot = spline_t(xi/xi(n), Omega_rot, deriv_type, dy_dx_a=0._WP, dy_dx_b=0._WP)

    ml%n_poly = n_poly
    ml%dt_Gamma_1 = Gamma_1
    ml%xi_1 = xi(n)
    ml%Omega_uni = 0._WP

    ml%uniform_rot = .FALSE.

    ! Finish

    return

  end function poly_model_t_

!****

  function V_2_1_ (this, x) result (V_2)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: V_2

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    ! Calculate V_2

    if (x /= 0._WP) then

       xi = x*this%xi_1

       Theta = this%sp_Theta%interp(xi)
       dTheta = this%sp_dTheta%interp(xi)

       V_2 = -(this%n_poly + 1._WP)*xi*dTheta/(Theta*x**2)

    else

       V_2 = (this%n_poly + 1._WP)*this%xi_1**2/3._WP

    endif

    ! Finish

    return

  end function V_2_1_

  
!****

  function V_2_v_ (this, x) result (V_2)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: V_2(SIZE(x))

    integer :: i

    ! Calculate V_2

    x_loop : do i = 1,SIZE(x)
       V_2(i) = this%V_2(x(i))
    end do x_loop

    ! Finish

    return

  end function V_2_v_

!****

  function As_1_ (this, x) result (As)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: As

    ! Calculate As

    As = this%V_2(x)*x**2 * &
         (this%n_poly/(this%n_poly + 1._WP) - 1._WP/this%Gamma_1(x))

    ! Finish

    return

  end function As_1_

!****

  function As_v_ (this, x) result (As)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: As(SIZE(x))

    integer :: i

    ! Calculate As

    x_loop : do i = 1,SIZE(x)
       As(i) = this%As(x(i))
    end do x_loop

    ! Finish

    return

  end function As_v_

!****

  function U_1_ (this, x) result (U)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: U

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    ! Calculate U

    xi = x*this%xi_1

    Theta = this%sp_Theta%interp(xi)
    dTheta = this%sp_dTheta%interp(xi)

    if(x /= 0._WP) then
       U = -xi*Theta**this%n_poly/dTheta
    else
       U = 3._WP
    endif

    ! Finish

    return

  end function U_1_

!****

  function U_v_ (this, x) result (U)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: U(SIZE(x))

    integer :: i

    ! Calculate U

    x_loop : do i = 1,SIZE(x)
       U(i) = this%U(x(i))
    end do x_loop

    ! Finish

    return

  end function U_v_

!****

  function D_1_ (this, x) result (D)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: D

    ! Calculate D = dlnrho/dlnx

    D = -this%V_2(x)*x**2*this%n_poly/(this%n_poly + 1._WP)
    
    ! Finish

    return

  end function D_1_

!****

  function D_v_ (this, x) result (D)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: D(SIZE(x))

    integer :: i

    ! Calculate D = dlnrho/dlnx

    x_loop : do i = 1,SIZE(x)
       D(i) = this%D(x(i))
    end do x_loop

    ! Finish

    return

  end function D_v_

!****

  function c_1_1_ (this, x) result (c_1)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: c_1

    real(WP) :: xi
    real(WP) :: dTheta
    real(WP) :: dTheta_1

    ! Calculate c_1

    xi = x*this%xi_1

    dTheta = this%sp_dTheta%interp(xi)
    dTheta_1 = this%sp_dTheta%interp(this%xi_1)

    if(x /= 0._WP) then
       c_1 = x*dTheta_1/dTheta
    else
       c_1 = -3._WP*dTheta_1/this%xi_1
    endif

    ! Finish

    return

  end function c_1_1_

!****

  function c_1_v_ (this, x) result (c_1)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: c_1(SIZE(x))

    integer :: i

    ! Calculate c_1

    x_loop : do i = 1,SIZE(x)
       c_1(i) = this%c_1(x(i))
    end do x_loop

    ! Finish

    return
    
  end function c_1_v_

!****

  function Gamma_1_1_ (this, x) result (Gamma_1)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: Gamma_1

    ! Calculate Gamma_1

    Gamma_1 = this%dt_Gamma_1

    ! Finish

    return

  end function Gamma_1_1_

!****
  
  function Gamma_1_v_ (this, x) result (Gamma_1)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: Gamma_1(SIZE(x))

    integer :: i

    ! Calculate Gamma_1
    
    x_loop : do i = 1,SIZE(x)
       Gamma_1(i) = this%Gamma_1(x(i))
    end do x_loop

    ! Finish

    return

  end function Gamma_1_v_

!****

  function nabla_ad_1_ (this, x) result (nabla_ad)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: nabla_ad

    ! Calculate nabla_ad (assume ideal gas)

    nabla_ad = 2._WP/5._WP

    ! Finish

    return

  end function nabla_ad_1_

!****
  
  function nabla_ad_v_ (this, x) result (nabla_ad)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: nabla_ad(SIZE(x))

    integer :: i

    ! Calculate nabla_ad
    
    x_loop : do i = 1,SIZE(x)
       nabla_ad(i) = this%nabla_ad(x(i))
    end do x_loop

    ! Finish

    return

  end function nabla_ad_v_

!****

  function delta_1_ (this, x) result (delta)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: delta

    ! Calculate delta (assume ideal gas)

    delta = 1._WP

    ! Finish

    return

  end function delta_1_

!****
  
  function delta_v_ (this, x) result (delta)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: delta(SIZE(x))

    integer :: i

    ! Calculate delta
    
    x_loop : do i = 1,SIZE(x)
       delta(i) = this%delta(x(i))
    end do x_loop

    ! Finish

    return

  end function delta_v_

!****

  $define $PROC $sub

  $local $NAME $1

  function ${NAME}_1_ (this, x) result ($NAME)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: $NAME

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_1_

!****

  function ${NAME}_v_ (this, x) result ($NAME)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: $NAME(SIZE(x))

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_v_

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
 
  function Omega_rot_1_ (this, x) result (Omega_rot)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    real(WP)                        :: Omega_rot

    ! Interpolate Omega_rot. If uniform_rot is .TRUE., use the uniform
    ! rate given by Omega_uni instead

    if (this%uniform_rot) then
       Omega_rot = this%Omega_uni
    else
       Omega_rot = this%sp_Omega_rot%interp(x)
    endif

    ! Finish

    return

  end function Omega_rot_1_

!****

  function Omega_rot_v_ (this, x) result (Omega_rot)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: Omega_rot(SIZE(x))

    ! Interpolate Omega_rot. If uniform_rot is .TRUE., use the uniform
    ! rate given by Omega_uni instead

    if (this%uniform_rot) then
       Omega_rot = this%Omega_uni
    else
       Omega_rot = this%sp_Omega_rot%interp(x)
    endif

    ! Finish

    return

  end function Omega_rot_v_

!****

  function pi_c_ (this) result (pi_c)

    class(poly_model_t), intent(in) :: this
    real(WP)                        :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c =  (this%n_poly + 1._WP)*this%xi_1**2/3._WP

    ! Finish

    return

  end function pi_c_

!****

  function is_zero_ (this, x) result (is_zero)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x
    logical                         :: is_zero

    real(WP) :: xi

    ! Determine whether the point at x has a vanishing pressure and/or density

    xi = x*this%xi_1

    is_zero = this%sp_Theta%interp(xi) == 0._WP

    ! Finish

    return

  end function is_zero_

!****

  subroutine attach_cache_ (this, cc)

    class(poly_model_t), intent(inout)    :: this
    class(cocache_t), pointer, intent(in) :: cc

    ! Enable a coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine attach_cache_

!****

  subroutine detach_cache_ (this)

    class(poly_model_t), intent(inout) :: this

    ! Detach the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine detach_cache_

!****

  subroutine fill_cache_ (this, x)

    class(poly_model_t), intent(inout) :: this
    real(WP), intent(in)               :: x(:)

    ! Fill the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine fill_cache_

!****

  $if ($MPI)

  subroutine bcast_ (ml, root_rank)

    type(poly_model_t), intent(inout) :: ml
    integer, intent(in)               :: root_rank

    ! Broadcast the poly_model_t

    call bcast(ml%sp_Theta, root_rank)
    call bcast(ml%sp_dTheta, root_rank)

    call bcast(ml%n_poly, root_rank)
    call bcast(ml%dt_Gamma_1, root_rank)
    call bcast(ml%xi_1, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

end module gyre_poly_model
