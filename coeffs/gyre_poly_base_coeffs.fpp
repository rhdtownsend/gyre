! Module   : gyre_poly_base_coeffs
! Purpose  : base structure coefficients for polytropic models
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

module gyre_poly_base_coeffs

  ! Uses

  use core_kinds
  use core_parallel
  use core_spline

  use gyre_base_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1
    procedure :: ${NAME}_v
  $endsub

  type, extends(base_coeffs_t) :: poly_base_coeffs_t
     private
     type(spline_t)   :: sp_Theta
     type(spline_t)   :: sp_dTheta
     real(WP)         :: dt_Gamma_1
     real(WP), public :: n_poly
     real(WP), public :: xi_1
   contains
     private
     procedure, public :: init
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     procedure, public :: pi_c
  end type poly_base_coeffs_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_bc
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: poly_base_coeffs_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, xi, Theta, dTheta, n_poly, Gamma_1, deriv_type)

    class(poly_base_coeffs_t), intent(out) :: this
    real(WP), intent(in)                   :: xi(:)
    real(WP), intent(in)                   :: Theta(:)
    real(WP), intent(in)                   :: dTheta(:)
    real(WP), intent(in)                   :: n_poly
    real(WP), intent(in)                   :: Gamma_1
    character(LEN=*), intent(in)           :: deriv_type

    integer  :: n
    real(WP) :: d2Theta(SIZE(xi))

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))

    ! Initialize the base_coeffs

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

    call this%sp_Theta%init(xi, Theta, dTheta)
    call this%sp_dTheta%init(xi, dTheta, d2Theta)

    this%n_poly = n_poly
    this%dt_Gamma_1 = Gamma_1
    this%xi_1 = xi(n)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bc (bc, root_rank)

    class(poly_base_coeffs_t), intent(inout) :: bc
    integer, intent(in)                      :: root_rank

    ! Broadcast the base_coeffs

    call bcast(bc%sp_Theta, root_rank)
    call bcast(bc%sp_dTheta, root_rank)

    call bcast(bc%n_poly, root_rank)
    call bcast(bc%dt_Gamma_1, root_rank)
    call bcast(bc%xi_1, root_rank)

    ! Finish

    return

  end subroutine bcast_bc

  $endif

!****

  function V_1 (this, x) result (V)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: V

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    ! Calculate V

    xi = x*this%xi_1

    Theta = this%sp_Theta%interp(xi)
    dTheta = this%sp_dTheta%interp(xi)

    V = -(this%n_poly + 1._WP)*xi*dTheta/Theta

    ! Finish

    return

  end function V_1

!****

  function V_v (this, x) result (V)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: V(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: As

    ! Calculate As

    As = this%V(x)*(this%n_poly/(this%n_poly + 1._WP) - 1._WP/this%Gamma_1(x))

    ! Finish

    return

  end function As_1

!****

  function As_v (this, x) result (As)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: As(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: U

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

  end function U_1

!****

  function U_v (this, x) result (U)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: U(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: c_1

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

  end function c_1_1

!****

  function c_1_v (this, x) result (c_1)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: c_1(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: Gamma_1

    ! Calculate Gamma_1

    Gamma_1 = this%dt_Gamma_1

    ! Finish

    return

  end function Gamma_1_1

!****
  
  function Gamma_1_v (this, x) result (Gamma_1)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: Gamma_1(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: nabla_ad

    ! Calculate nabla_ad (assume ideal gas)

    nabla_ad = 2._WP/5._WP

    ! Finish

    return

  end function nabla_ad_1

!****
  
  function nabla_ad_v (this, x) result (nabla_ad)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: nabla_ad(SIZE(x))

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

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: delta

    ! Calculate delta (assume ideal gas)

    delta = 1._WP

    ! Finish

    return

  end function delta_1

!****
  
  function delta_v (this, x) result (delta)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: delta(SIZE(x))

    integer :: i

    ! Calculate delta
    
    x_loop : do i = 1,SIZE(x)
       delta(i) = this%delta(x(i))
    end do x_loop

    ! Finish

    return

  end function delta_v

!****

  function pi_c (this)

    class(poly_base_coeffs_t), intent(in) :: this
    real(WP)                              :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c =  (this%n_poly + 1._WP)*this%xi_1**2/3._WP

    ! Finish

    return

  end function pi_c

end module gyre_poly_base_coeffs
