! Module   : gyre_poly_seg
! Purpose  : polytropic model segment
!
! Copyright 2016 Rich Townsend
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

module gyre_poly_seg

  ! Uses

  use core_kinds

  use gyre_spline
  use gyre_model_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure, public :: ${NAME}
  $endsub

  type :: poly_seg_t
     private
     type(r_spline_t) :: sp_Theta
     type(r_spline_t) :: sp_dTheta
     real(WP)         :: mu_i
     real(WP)         :: mu_s
     real(WP)         :: B
     real(WP)         :: v_i
     real(WP)         :: xi_s
     real(WP)         :: n_poly
     real(WP)         :: Gamma_1_
   contains
     private
     procedure :: mu
     $PROC_DECL(V_2)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(dU)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(delta)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(dnabla_ad)
     $PROC_DECL(Omega_rot)
     $PROC_DECL(dOmega_rot)
  end type poly_seg_t
 
  ! Interfaces

  interface poly_seg_t
     module procedure poly_seg_t_
  end interface poly_seg_t

  ! Access specifiers

  private

  public :: poly_seg_t

  ! Procedures

contains

  function poly_seg_t_ (x, Theta, dTheta, mu_i, mu_s, B, xi_s, n_poly, Gamma_1) result (ps)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: Theta(:)
    real(WP), intent(in) :: dTheta(:)
    real(WP), intent(in) :: mu_i
    real(WP), intent(in) :: mu_s
    real(WP), intent(in) :: B
    real(WP), intent(in) :: xi_s
    real(WP), intent(in) :: n_poly
    real(WP), intent(in) :: Gamma_1
    type(poly_seg_t)     :: ps

    real(WP) :: xi(SIZE(x))
    real(WP) :: d2Theta(SIZE(x))

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    ! Construct the poly_seg_t

    xi = x*xi_s

    if (n_poly /= 0._WP) then

       where (xi /= 0._WP)
          d2Theta = -2._WP*dTheta/xi - B*Theta**n_poly
       elsewhere
          d2Theta = -1._WP/3._WP
       end where

    else

       d2Theta = -1._WP/3._WP

    endif

    ps%sp_Theta = r_spline_t(x, Theta, dTheta*xi_s)
    ps%sp_dTheta = r_spline_t(x, dTheta, d2Theta*xi_s)

    ps%mu_i = mu_i
    ps%mu_s = mu_s

    ps%B = B
    ps%v_i = xi(1)**2*dTheta(1)

    ps%xi_s = xi_s

    ps%n_poly = n_poly
    ps%Gamma_1_ = Gamma_1

    ! Finish

    return

  end function poly_seg_t_

  !****

  function mu (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: mu

    real(WP) :: xi
    real(WP) :: v

    ! Calculate the mass coordinate mu

    xi = x*this%xi_s
    v = xi**2*this%sp_dTheta%f(x)

    mu = this%mu_i - (v - this%v_i)/this%B

    ! Finish

    return

  end function mu

  !****

  function V_2 (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: V_2

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    ! Calculate V_2

    if (x /= 0._WP) then

       xi = x*this%xi_s

       Theta = this%sp_Theta%f(x)
       dTheta = this%sp_dTheta%f(x)

       if (Theta == 0._WP) Theta = TINY(0._WP)

       V_2 = -(this%n_poly + 1._WP)*xi*dTheta/(Theta*x**2)

    else

       V_2 = (this%n_poly + 1._WP)*this%xi_s**2/3._WP

    endif

    ! Finish

    return

  end function V_2

  !****

  function As (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: As

    ! Calculate As

    As = this%V_2(x)*x**2 * &
         (this%n_poly/(this%n_poly + 1._WP) - 1._WP/this%Gamma_1_)

    ! Finish

    return

  end function As

  !****

  function U (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: U

    real(WP) :: xi
    real(WP) :: Theta

    ! Calculate U

    if (x /= 0._WP) then

       xi = x*this%xi_s

       Theta = this%sp_Theta%f(x)

       U = xi**3*Theta**this%n_poly/this%mu(x)

    else

       U = 3._WP

    endif

    ! Finish

    return

  end function U

  !****

  function dU (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: dU

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    ! Calculate dU

    if (x /= 0._WP) then

       xi = x*this%xi_s

       Theta = this%sp_Theta%f(x)
       dTheta = this%sp_dTheta%f(x)
       
       if (Theta == 0._WP) Theta = TINY(0._WP)

       dU = 3._WP + this%n_poly*xi*dTheta/Theta - xi**3*Theta**this%n_poly/this%mu(x)

    else

       dU = 0._WP

    endif

    ! Finish

    return

  end function dU

  !****

  function c_1 (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: c_1

    ! Calculate c_1

    if (x /= 0._WP) then
       c_1 = x**3/(this%mu(x)/this%mu_s)
    else
       c_1 = 3._WP*this%mu_s/this%xi_s**3
    endif

    ! Finish

    return

  end function c_1

  !****

  function Gamma_1 (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: Gamma_1

    ! Calculate Gamma_1

    Gamma_1 = this%Gamma_1_

    ! Finish

    return

  end function Gamma_1

  !****

  function delta (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: delta

    ! Calculate delta (assume ideal gas)

    delta = 1._WP

    ! Finish

    return

  end function delta

  !****

  function nabla_ad (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: nabla_ad

    ! Calculate nabla_ad (assume ideal gas)

    nabla_ad = 2._WP/5._WP

    ! Finish

    return

  end function nabla_ad

  !****

  function dnabla_ad (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: dnabla_ad

    ! Calculate dlnnabla_ad/dlnx (assume ideal gas)

    dnabla_ad = 0._WP

    ! Finish

    return

  end function dnabla_ad

  !****

  function Omega_rot (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: Omega_rot

    ! Calculate Omega_rot

    Omega_rot = 0._WP

    ! Finish

    return

  end function Omega_rot
  
  !****

  function dOmega_rot (this, x)

    class(poly_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: dOmega_rot

    ! Calculate dlnOmega_rot/dlnx

    dOmega_rot = 0._WP

    ! Finish

    return

  end function dOmega_rot
  
end module gyre_poly_seg
