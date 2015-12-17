! Module   : gyre_evol_model
! Purpose  : stellar model (evolutionary)
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_evol_model

  ! Uses

  use core_kinds
  
  use gyre_constants
  use gyre_coords
  use gyre_evol_model_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  $define $SET_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME} => set_${NAME}_
  $endsub

  type, extends (model_t) :: evol_model_t
     private
     type(coords_t), allocatable    :: co(:)
     type(evol_seg_t), allocatable  :: es(:)
     real(WP), public               :: M_star
     real(WP), public               :: R_star
     real(WP), public               :: L_star
     integer                        :: n_co
   contains
     private
     $SET_DECL(V_2)
     $SET_DECL(As)
     $SET_DECL(U)
     $SET_DECL(c_1)
     $SET_DECL(Gamma_1)
     $SET_DECL(delta)
     $SET_DECL(nabla_ad)
     $SET_DECL(nabla)
     $SET_DECL(beta_rad)
     $SET_DECL(c_rad)
     $SET_DECL(c_thm)
     $SET_DECL(c_dif)
     $SET_DECL(c_eps_ad)
     $SET_DECL(c_eps_S)
     $SET_DECL(kappa_ad)
     $SET_DECL(kappa_S)
     $SET_DECL(Omega_rot)
     $PROC_DECL(V_2)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(dU)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(delta)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(dnabla_ad)
     $PROC_DECL(nabla)
     $PROC_DECL(beta_rad)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(Omega_rot)
     $PROC_DECL(dOmega_rot)
     $PROC_DECL(M_r)
     $PROC_DECL(P)
     $PROC_DECL(rho)
     $PROC_DECL(T)
     procedure, public :: scaffold => scaffold_
     procedure, public :: delta_p => delta_p_
     procedure, public :: delta_g => delta_g_
  end type evol_model_t
 
  ! Interfaces

  interface evol_model_t
     module procedure evol_model_t_
  end interface evol_model_t

  ! Access specifiers

  private

  public :: evol_model_t

  ! Procedures

contains

  function evol_model_t_ (ml_p, M_star, R_star, L_star, x) result (ml)

    type(model_par_t), intent(in) :: ml_p
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    real(WP), intent(in)          :: x(:)
    type(evol_model_t)            :: ml

    integer :: s

    ! Construct the evol_model_t

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ml%n_co = SIZE(x)

    allocate(ml%co(ml%n_co))
    call partition(x, ml%co, ml%n_s)

    allocate(ml%es(ml%n_s))

    seg_loop : do s = 1, ml%n_s
       ml%es(s) = evol_seg_t(ml_p, PACK(x, MASK=ml%co%s == s))
    end do seg_loop

    ! Finish

    return

  end function evol_model_t_

  !****

  $define $SET $sub

  $local $NAME $1

  subroutine set_${NAME}_ (this, f)

    class(evol_model_t), intent(inout) :: this
    real(WP), intent(in)               :: f(:)

    integer :: s

    $CHECK_BOUNDS(SIZE(f),this%n_co)

    ! Set the data for $NAME

    seg_loop : do s = 1, this%n_s
       call this%es(s)%set_$NAME(PACK(f, MASK=this%co%s == s))
    end do seg_loop

    ! Finish

    return

  end function set_${NAME}_

  $endif

  $SET(V_2)
  $SET(A_s)
  $SET(U)
  $SET(c_1)
  $SET(Gamma_1)
  $SET(delta)
  $SET(nabla_ad)
  $SET(nabla)
  $SET(beta_rad)
  $SET(c_rad)
  $SET(c_thm)
  $SET(c_dif)
  $SET(c_eps_ad)
  $SET(c_eps_S)
  $SET(kappa_ad)
  $SET(kappa_S)
  $SET(Omega_rot)

  !****

  $define $PROC_1 $sub

  $local $NAME $1

  function ${NAME}_1_ (this, co) result ($NAME)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: $NAME

    ! Evaluate $NAME

    associate (s => co%s, x => co%x)

      $NAME = this%es(s)%$NAME(x)

    end associate

    ! Finish

    return

  end function ${NAME}_1_

  $endsub

  $PROC_1(V_2)
  $PROC_1(As)
  $PROC_1(U)
  $PROC_1(dU)
  $PROC_1(c_1)
  $PROC_1(Gamma_1)
  $PROC_1(delta)
  $PROC_1(nabla_ad)
  $PROC_1(dnabla_ad)
  $PROC_1(nabla)
  $PROC_1(beta_rad)
  $PROC_1(c_rad)
  $PROC_1(dc_rad)
  $PROC_1(c_thm)
  $PROC_1(c_dif)
  $PROC_1(c_eps_ad)
  $PROC_1(c_eps_S)
  $PROC_1(kappa_ad)
  $PROC_1(kappa_S)
  $PROC_1(Omega_rot)
  $PROC_1(dOmega_rot)

  !****

  function M_r_1_ (this, mc) result (M_r)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: M_r

    ! Evaluate the fractional mass coordinate

    M_r = this%M_star*(co%x**3*/this%c_1(co))

    ! Finish

    return

  end function M_r_1_
    
  !****

  function P_1_ (this, co) result (P)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: P

    ! Evaluate the total pressure

    P = (G_GRAVITY*this%M_star/(4._WP*PI*this%R_star**4))*&
        (co%x**2*this%U(co)/(this%c_1(co)*this%V(co)))

    ! Finish

    return

  end function P_1_
    
  !****

  function rho_1_ (this, co) result (rho)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: rho

    ! Evaluate the density

    rho = (this%M_star/(4._WP*PI*this%R_star)**3*(this%U(co)/this%c_1(co))

    ! Finish

    return

  end function rho_1_
    
  !****

  function T_1_ (this, co) result (T)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: T

    ! Evaluate the temperature

    T = (3._WP*this%beta_r(co)*this%P(co)/A_RADIATION)**0.25_WP

    ! Finish

    return

  end function T_1_
    
  !****

  $define $PROC_V $sub

  $local $NAME $1

  function ${NAME}_v_ (this, co) result($NAME)

    class(evol_model_t), intent(in) :: this
    type(coords_t), intent(in)      :: co
    real(WP)                        :: $NAME(SIZE(co))

    integer :: i

    ! Evaluate $NAME

    !$OMP PARALLEL DO
    coords_loop : do i = 1, SIZE(co)
       $NAME(i) = this%$NAME(co(i))
    end do coords_loop

    ! Finish

    return

  end function ${NAME}_v_

  $endsub

  $PROC_V(V_2)
  $PROC_V(As)
  $PROC_V(U)
  $PROC_V(dU)
  $PROC_V(c_1)
  $PROC_V(Gamma_1)
  $PROC_V(delta)
  $PROC_V(nabla_ad)
  $PROC_V(dnabla_ad)
  $PROC_V(nabla)
  $PROC_V(beta_rad)
  $PROC_V(c_rad)
  $PROC_V(dc_rad)
  $PROC_V(c_thm)
  $PROC_V(c_dif)
  $PROC_V(c_eps_ad)
  $PROC_V(c_eps_S)
  $PROC_V(kappa_ad)
  $PROC_V(kappa_S)
  $PROC_V(Omega_rot)
  $PROC_V(dOmega_rot)
  $PROC_V(M_r)
  $PROC_V(P)
  $PROC_V(rho)
  $PROC_V(T)

  !****

  function scaffold_ (this) result (co)

    class(evol_model_t), intent(in) :: this
    class(coords_t), allocatable    :: co(:)

    ! Return the grid scaffold

    co = this%co

    ! Finish

    return

  end function scaffold_

  !****

  function delta_p_ (this) result (delta_p)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_p

    real(WP) :: V_2(this%n_co)
    real(WP) :: c_1(this%n_co)
    real(WP) :: Gamma_1(this%n_co)
    real(WP) :: f(this%n_co)

    ! Calculate the p-mode (large) frequency separation

    V_2 = this%V_2(this%co)
    c_1 = this%c_1(this%co)
    Gamma_1 = this%Gamma_1(this%co)

    f = Gamma_1/(c_1*V_2)

    delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/ &
              integrate(this%co%x, f)
       
    ! Finish

    return

  end function delta_p_

  !****

  function delta_g_(this) result(delta_g)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_g

    real(WP) :: As(this%n_co)
    real(WP) :: c_1(this%n_co)
    real(WP) :: f(this%n_co)

    ! Calculate the g-mode inverse period separation

    As = this%As(this%co)
    c_1 = this%c_1(this%co)

    where (this%co%x /= 0._WP)
       f = MAX(As/c_1, 0._WP))/this%co%x
    elsewhere
       f = 0._WP
    end where

    delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/PI**2* &
              integrate(this%co%x, f)

    ! Finish

    return

  end function delta_g_

end module gyre_evol_model
