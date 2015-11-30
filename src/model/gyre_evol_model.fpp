! Module   : gyre_evol_model
! Purpose  : stellar model (evolutionary)
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

module gyre_evol_model

  ! Uses

  use core_kinds
  
  use gyre_constants
  use gyre_evol_seg
  use gyre_model
  use gyre_model_par
  use gyre_part

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME} => set_${NAME}_
  $endsub

  type, extends(model_t) :: evol_model_t
     private
     type(model_par_t)             :: ml_p
     type(evol_seg_t), allocatable :: sg(:)
     real(WP), allocatable         :: x(:)
     real(WP), public              :: M_star
     real(WP), public              :: R_star
     real(WP), public              :: L_star
   contains
     private
     procedure, public :: seg => seg_
     procedure, public :: n_seg => n_seg_
     $PROC_DECL(V_2)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(delta)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(nabla)
     $PROC_DECL(c_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(Omega_rot)
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

  function evol_model_t_(ml_p, M_star, R_star, L_star, x)

    type(model_par_t), intent(in) :: ml_p
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    real(WP), intent(in)          :: x(:)

    real(WP), allocatable :: x_c(:)
    type(part_t)          :: pt
    integer               :: s

    ! Construct the evol_model_t

    ml%ml_p = ml_p

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ! If necessary, add the central point to x

    if(ml_p%add_center .AND. x(1) /= 0._WP) then
       ml%x = [0._WP,x]
    else
       ml%x = x
    endif

    ! Create segments

    pt = part_t(ml%x)

    seg_loop : do s = 1, pt%n_seg()
       ml%sg(s) = evol_seg_t(ml%x(pt%i_seg(s)))
    end do seg_loop

    end associate

    ! Finish

    return

  end function evol_model_t_

  !****

  function seg_(this, s) result(sg)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    class(seg_t), pointer           :: sg

    ! Return a pointer to the s'th segment

    sg => this%sg(s)

    ! Finish

    return

  end function seg_

  !****

  function n_seg_(this) result(n_seg)

    class(evol_model_t), intentt(in) :: this
    integer                          :: n_seg

    ! Return the number of segments

    n_seg = SIZE(this%sg)

    ! Finish

  end function n_seg_

  !****

  $define $SET_CONT $sub

  $local $NAME $1

  function set_${NAME}_(this, y) result(sp)

    class(evol_model_t), intent(in) :: this
    real(WP), intent(in)            :: y(:)

    type(part_t)          :: pt
    integer, allocatable  :: i_dbl(:)
    integer               :: s
    integer, allocatable  :: i_unq(:)
    real(WP), allocatable :: x_unq(:)
    real(WP), allocatable :: y_unq(:)
    real(WP), allocatable :: dy_dx_unq(:)
    real(WP), allocatable :: dy_dx(:)
    integer, allocatable  :: i_seg(:)

    $CHECK_BOUNDS(y,this%x)

    ! Set up data for the continuous function $NAME in the s'th
    ! segment

    pt = part_t(this%x)

    ! Verify continuity

    i_dbl = pt%i_dbl()

    cont_loop : do s = 1, SIZE(i_dbl)
       $ASSERT(y(i_dbl(s)) == y(i_dbl(s)+1),Discontinuous $NAME at double point)
    end do cont_loop

    ! Calculate derivatives globally

    i_unq = pt%i_unq()

    x_unq = x(i_unq)
    y_unq = y(i_unq)

    if(x_unq(1) == 0._WP) then
       dy_dx_unq = eval_dy_dx(x_unq, y_unq, dy_dx_a=0._WP)
    else
       dy_dx_unq = eval_dy_dx(x_unq, y_unq)
    endif

    dy_dx(i_unq) = dy_dx_unq
    dy_dx(i_dbl+1) = dy_dx(i_dbl)

    ! Loop through the segments

    seg_loop : do s = 1, this%n_seg()

       i_seg = pt%i_seg(s)

       ! Store segment data

       this%sg(s)%set_${NAME}(interp_t(x(i_seg), y(i_seg), dy_dx(i_seg)))

    end do seg_loop

    ! Finish

    return

  end function set_${NAME}

  $endif

  $SET_CONT(c_1)

  !****

  $define $SET_DISCONT $sub

  $local $NAME $1

  function set_${NAME}_(this, y) result(sp)

    class(evol_model_t), intent(in) :: this
    real(WP), intent(in)            :: y(:)

    type(part_t)          :: pt
    integer               :: s
    integer, allocatable  :: i_seg(:)
    real(WP), allocatable :: x_seg(:)
    real(WP), allocatable :: y_seg(:)
    real(WP), allocatable :: dy_dx_seg(:)

    $CHECK_BOUNDS(SIZE(y),SIZE(this%x))

    ! Set up data for the discontinuous function $NAME in the s'th
    ! segment

    pt = part_t(this%x)

    ! Loop through the segments

    seg_loop : do s = 1, this%n_seg()

       i_seg = pt%i_seg(s)

       ! Calculate derivatives locally

       x_seg = x(i_seg)
       y_seg = y(i_seg)

       if(x_seg(1) == 0._WP) then
          dy_dx_seg = eval_dy_dx(x_seg, y_seg, dy_dx_a=0._WP)
       else
          dy_dx_seg = eval_dy_dx(x_seg, y_seg)
       endif

       ! Store segment data
       
       this%sg(s)%set_${NAME}(interp_t(x_seg, y_seg, dy_dx_seg))

    end do seg_loop

    ! Finish

    return

  end function set_${NAME}_

  $endif

  $SET_DISCONT(V_2)
  $SET_DISCONT(A_s)
  $SET_DISCONT(U)
  $SET_DISCONT(Gamma_1)
  $SET_DISCONT(delta)
  $SET_DISCONT(nabla_ad)
  $SET_DISCONT(nabla)
  $SET_DISCONT(c_rad)
  $SET_DISCONT(c_thm)
  $SET_DISCONT(c_dif)
  $SET_DISCONT(c_eps_ad)
  $SET_DISCONT(c_eps_S)
  $SET_DISCONT(kappa_ad)
  $SET_DISCONT(kappa_S)
  $SET_DISCONT(Omega_rot)

  !****

  function delta_p_(this) result(delta_p)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_p

    real(WP)                  :: int_f
    integer                   :: s
    type(evol_seg_t), pointer :: sg
    real(WP), allocatable     :: V_2(:)
    real(WP), allocatable     :: c_1(:)
    real(WP), allocatable     :: Gamma_1(:)
    real(WP), allocatable     :: f(:)

    ! Calculate the p-mode (large) frequency separation

    int_f = 0._WP

    seg_loop : do s = 1, this%n_seg()

       sg => this%seg(s)

       V_2 = sg%V_2(sg%x)
       c_1 = sg%c_1(sg%x)
       Gamma_1 = sg%Gamma_1(sg%x)

       f = Gamma_1/(c_1*V_2)

       int_f = int_f + integrate(sg%x, f)

    end do seg_loop

    delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/int_f
       
    ! Finish

    return

  end function delta_p_

  !****

  function delta_g_(this) result(delta_g)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_g

    real(WP)                  :: int_f
    integer                   :: s
    type(evol_seg_t), pointer :: sg
    real(WP), allocatable     :: As(:)
    real(WP), allocatable     :: c_1(:)
    real(WP), allocatable     :: f(:)

    ! Calculate the g-mode inverse period separation

    int_f = 0._WP

    seg_loop : do s = 1, this%n_seg()

       sg => this%seg(s)

       As = sg%As(sg%x)
       c_1 = sg%c_1(sg%x)

       where(sg%x /= 0._WP)
          f = MAX(As/c_1, 0._WP))/sg%x
       elsewhere
          f = 0._WP
       end where

       int_f = int_f + integrate(sg%x, f)

    end do seg_loop

    delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)*int_f/PI**2

    ! Finish

    return

  end function delta_g_

end module gyre_evol_model
