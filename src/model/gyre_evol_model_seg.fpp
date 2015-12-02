! Module   : gyre_evol_model_seg
! Purpose  : stellar model segment (evolutionary)
!
! Copyright 2015 Rich Townsend
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

module gyre_evol_model_seg

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_model_seg
  use gyre_interp
  use gyre_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $DATA_DECL $sub
    $local $NAME $1
    type(interp_t) :: ip_${NAME}
    logical        :: ${NAME}_def = .FALSE.
  $endif

  $define $SET_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME} => set_${NAME}_
  $endsub

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure         :: ${NAME}_1_
    procedure         :: ${NAME}_v_
  $endsub

  type, extends (model_seg_t) :: evol_model_seg_t
     private
     $DATA_DECL(V_2)
     $DATA_DECL(As)
     $DATA_DECL(U)
     $DATA_DECL(c_1)
     $DATA_DECL(Gamma_1)
     $DATA_DECL(delta)
     $DATA_DECL(nabla_ad)
     $DATA_DECL(nabla)
     $DATA_DECL(c_rad)
     $DATA_DECL(c_thm)
     $DATA_DECL(c_dif)
     $DATA_DECL(c_eps_ad)
     $DATA_DECL(c_eps_S)
     $DATA_DECL(kappa_ad)
     $DATA_DECL(kappa_S)
     $DATA_DECL(Omega_rot)
     $DATA_DECL(T)
     real(WP)                  :: M_star
     real(WP)                  :: R_star
     real(WP)                  :: L_star
     character(:), allocatable :: deriv_type
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
     $SET_DECL(c_rad)
     $SET_DECL(c_thm)
     $SET_DECL(c_dif)
     $SET_DECL(c_eps_ad)
     $SET_DECL(c_eps_S)
     $SET_DECL(kappa_ad)
     $SET_DECL(kappa_S)
     $SET_DECL(Omega_rot)
     $SET_DECL(T)
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
     $PROC_DECL(m)
     $PROC_DECL(P)
     $PROC_DECL(rho)
     $PROC_DECL(T)
  end type evol_model_seg_t
 
  ! Interfaces

  interface evol_model_seg_t
     module procedure evol_model_seg_t_part_
  end interface evol_model_seg_t

  ! Access specifiers

  private

  public :: evol_model_seg_t

  ! Procedures

contains

  function evol_model_seg_t_part_ (ml_p, M_star, R_star, L_star, x) result (ms)

    type(model_par_t), intent(in)       :: ml_p
    real(WP), intent(in)                :: M_star
    real(WP), intent(in)                :: R_star
    real(WP), intent(in)                :: L_star
    real(WP), intent(in)                :: x(:)
    type(evol_model_seg_t), allocatable :: ms(:)

    ! Construct an array of segments, bt splitting x at double points

    ms%seg_t = seg_t(x)

    ms%M_star = M_star
    ms%R_star = R_star
    ms%L_star = L_star

    seg_loop : do s = 1, SIZE(ms)
       ms(s)%deriv_type = ml_p%deriv_type
    end do seg_loop

    ! Finish

    return

  end function evol_model_seg_t_part_

  !****

  $def $SET $sub

  $local $NAME $1

  subroutine set_${NAME}_ (this, y)

    class(evol_model_seg_t), intent(inout) :: this
    real(WP), intent(in)                   :: y(:)

    $CHECK_BOUNDS(SIZE(y),this%n)

    ! Set the data for $NAME

    if (this%x(1) == 0._WP) then
       this%ip_${NAME} = interp_t(this%x, y, this%deriv_type, dy_dx_a=0._WP)
    else
       this%ip_${NAME} = interp_t(this%x, y, this%deriv_type)
    endif
    
    this%${NAME}_def = .TRUE.

    ! Finish

    return

  end subroutine set_${NAME}_

  $endif

  $SET(V_2)
  $SET(As)
  $SET(U)
  $SET(c_1)
  $SET(Gamma_1)
  $SET(delta)
  $SET(nabla_ad)
  $SET(nabla)
  $SET(c_rad)
  $SET(c_thm)
  $SET(c_dif)
  $SET(c_eps_ad)
  $SET(c_eps_S)
  $SET(kappa_ad)
  $SET(kappa_S)
  $SET(Omega_rot)

  !****

  $define $PROC $sub

  $local $NAME $1

  function ${NAME}_1_ (this, x) result ($NAME)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x
    real(WP)                            :: $NAME

    ! Interpolate $NAME

    if (this%${NAME}_def) then

       $NAME = this%ip_${NAME}%f(x)

    else

       $ABORT(No model data provided for $NAME)

    endif
       
    ! Finish

    return

  end function ${NAME}_1_

  function ${NAME}_v_(this, x) result($NAME)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x(:)
    real(WP)                            :: $NAME(SIZE(x))

    ! Interpolate $NAME

    if(this%${NAME}_def) then

       $NAME = this%ip_${NAME}%f(x)

    else

       $ABORT(No model data provided for $NAME)

    endif

    ! Finish

    return

  end function ${NAME}_v_

  $endsub

  $PROC(V_2)
  $PROC(As)
  $PROC(U)
  $PROC(c_1)
  $PROC(Gamma_1)
  $PROC(delta)
  $PROC(nabla_ad)
  $PROC(nabla)
  $PROC(c_rad)
  $PROC(c_thm)
  $PROC(c_dif)
  $PROC(c_eps_ad)
  $PROC(c_eps_S)
  $PROC(kappa_ad)
  $PROC(kappa_S)
  $PROC(Omega_rot)
  $PROC(T)

!****

  $define $DPROC $sub

  $local $NAME $1

  function d${NAME}_1_(this, x) result(d$NAME)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x
    real(WP)                            :: d$NAME

    if (this%${NAME}_def) then

       if(x > 0._WP) then
          d$NAME = x*this%ip_${NAME}%df_dx(x)/this%ip_${NAME}%f(x)
       else
          d$NAME = 0._WP
       endif

    else

       $ABORT(No model data provided for $NAME)

    endif

    ! Finish

    return

  end function d${NAME}_1_

  function d${NAME}_v_(this, x) result(d$NAME)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x(:)
    real(WP)                            :: d$NAME(SIZE(x))

    integer :: j

    ! Interpolate dln$NAME/dlnx

    if(this%${NAME}_def) then

       where(x > 0._WP)
          d$NAME = x*this%ip_${NAME}%df_dx(x)/this%ip_${NAME}%f(x)
       elsewhere
          d$NAME = 0._WP
       endwhere

    else

       $ABORT(No model data provided for $NAME)

    endif
       
    ! Finish

    return

  end function d${NAME}_v_

  $endsub

  $DPROC(nabla_ad)
  $DPROC(c_rad)
  $DPROC(Omega_rot)

  !****

  function M_r_1_(this, x) result(m)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x
    real(WP)                            :: M_r

    ! Calculate the fractional mass coordinate

    M_r = this%M_star*(x**3*/this%c_1(x))

    ! Finish

    return

  end function M_r_1_
    
  !****

  function M_r_v_(this, x) result(m)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x(:)
    real(WP)                            :: M_r(SIZE(x))

    ! Calculate the fractional mass coordinate

    M_r = this%M_star*(x**3*/this%c_1(x))

    ! Finish

    return

  end function M_r_v_
    
  !****

  function P_1_(this, x) result(P)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x
    real(WP)                            :: P

    ! Calculate the pressure

    P = (G_GRAVITY*this%M_star/(4._WP*PI*this%R_star**4))*&
        (x**2*this%U(x)/(this%c_1(x)*this%V(x)))

    ! Finish

    return

  end function P_1_
    
  !****

  function P_v_(this, x) result(P)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x(:)
    real(WP)                            :: P(SIZE(x))

    ! Calculate the pressure

    P = (G_GRAVITY*this%M_star/(4._WP*PI*this%R_star**4))*&
        (x**2*this%U(x)/(this%c_1(x)*this%V(x)))

    ! Finish

    return

  end function P_v_
    
  !****

  function rho_1_(this, x) result(P)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x
    real(WP)                            :: P

    ! Calculate the density

    rho = (this%M_star/(4._WP*PI*this%R_star)**3*(this%U(x)/this%c_1(x))

    ! Finish

    return

  end function rho_1_
    
  !****

  function rho_v_(this, x) result(rho)

    class(evol_model_seg_t), intent(in) :: this
    real(WP), intent(in)                :: x(:)
    real(WP)                            :: rho(SIZE(x))

    ! Calculate the density

    rho = (this%M_star/(4._WP*PI*this%R_star)**3*(this%U(x)/this%c_1(x))

    ! Finish

    return

  end function rho_v_
    
end module gyre_evol_model_seg
