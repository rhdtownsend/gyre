! Module   : gyre_evol_seg
! Purpose  : evolutionary model segment
!
! Copyright 2015-2016 Rich Townsend
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

module gyre_evol_seg

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_spline
  use gyre_model_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $DATA_DECL $sub
    $local $NAME $1
    type(r_spline_t) :: sp_${NAME}
    logical          :: df_${NAME} = .FALSE.
  $endsub

  $define $SET_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME}
  $endsub

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure, public :: ${NAME}
  $endsub

  type :: evol_seg_t
     private
     $DATA_DECL(V_2)
     $DATA_DECL(As)
     $DATA_DECL(U)
     $DATA_DECL(c_1)
     $DATA_DECL(Gamma_1)
     $DATA_DECL(delta)
     $DATA_DECL(nabla_ad)
     $DATA_DECL(nabla)
     $DATA_DECL(beta_rad)
     $DATA_DECL(c_rad)
     $DATA_DECL(c_thm)
     $DATA_DECL(c_dif)
     $DATA_DECL(c_eps_ad)
     $DATA_DECL(c_eps_S)
     $DATA_DECL(kappa_ad)
     $DATA_DECL(kappa_S)
     $DATA_DECL(Omega_rot)
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
  end type evol_seg_t

  
  ! Interfaces

  interface evol_seg_t
     module procedure evol_seg_t_
  end interface evol_seg_t

  ! Access specifiers

  private

  public :: evol_seg_t

  ! Procedures

contains

  function evol_seg_t_ (ml_p) result (es)

    type(model_par_t), intent(in) :: ml_p
    type(evol_seg_t)              :: es

    ! Construct the evol_seg_t

    es%deriv_type = ml_p%deriv_type

    ! Finish

    return

  end function evol_seg_t_

  !****

  $define $SET $sub

  $local $NAME $1

  subroutine set_${NAME} (this, x, f)

    class(evol_seg_t), intent(inout) :: this
    real(WP), intent(in)             :: x(:)
    real(WP), intent(in)             :: f(:)

    $CHECK_BOUNDS(SIZE(f),SIZE(x))

    ! Set the data for $NAME

    if (x(1) == 0._WP) then
       this%sp_${NAME} = r_spline_t(x, f, this%deriv_type, df_dx_a=0._WP)
    else
       this%sp_${NAME} = r_spline_t(x, f, this%deriv_type)
    endif
    
    this%df_${NAME} = .TRUE.

    ! Finish

    return

  end subroutine set_${NAME}

  $endsub

  $SET(V_2)
  $SET(As)
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

  $define $PROC $sub

  $local $NAME $1

  function ${NAME} (this, x)

    class(evol_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: ${NAME}

    ! Interpolate $NAME

    if (this%df_${NAME}) then

       $NAME = this%sp_${NAME}%f(x)

    else

       $ABORT(No model data provided for $NAME)

    endif
       
    ! Finish

    return

  end function ${NAME}
  
  $endsub

  $PROC(V_2)
  $PROC(As)
  $PROC(U)
  $PROC(c_1)
  $PROC(Gamma_1)
  $PROC(delta)
  $PROC(nabla_ad)
  $PROC(nabla)
  $PROC(beta_rad)
  $PROC(c_rad)
  $PROC(c_thm)
  $PROC(c_dif)
  $PROC(c_eps_ad)
  $PROC(c_eps_S)
  $PROC(kappa_ad)
  $PROC(kappa_S)
  $PROC(Omega_rot)

  !****

  $define $DPROC $sub

  $local $NAME $1

  function d${NAME} (this, x)

    class(evol_seg_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP)                      :: d${NAME}

    ! Interpolate dln${NAME}/dlnx

    if (this%df_${NAME}) then

       if (x > 0._WP) then
          d$NAME = x*this%sp_${NAME}%df_dx(x)/this%sp_${NAME}%f(x)
       else
          d$NAME = 0._WP
       endif

    else

       $ABORT(No model data provided for $NAME)

    endif

    ! Finish

    return

  end function d${NAME}

  $endsub

  $DPROC(U)
  $DPROC(nabla_ad)
  $DPROC(c_rad)
  $DPROC(Omega_rot)

end module gyre_evol_seg
