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
  use gyre_evol_model_seg
  use gyre_grid_part
  use gyre_model_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME} => set_${NAME}_
  $endsub

  type, extends (model_t) :: evol_model_t
     private
     real(WP), public :: M_star
     real(WP), public :: R_star
     real(WP), public :: L_star
   contains
     private
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
     $PROC_DECL(T)
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

    integer      :: s

    ! Construct the evol_model_t

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ml%ms = evol_model_seg_t(ml_p, M_star, R_star, L_star, x)

    ml%n_s = SIZE(ml%ms)

    ! Finish

    return

  end function evol_model_t_

  !****

  $define $SET $sub

  $local $NAME $1

  subroutine set_${NAME}_ (this, y)

    class(evol_model_t), intent(inout) :: this
    real(WP), intent(in)               :: y(:)

    integer :: s

    ! Set the data for $NAME

    seg_loop : do s = 1, this%n_s

       associate (ms => this%ms(s))
         call ms%set_${NAME}(y(ms%i_a:ms%i_b))
       end associate

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
  $SET(c_rad)
  $SET(c_thm)
  $SET(c_dif)
  $SET(c_eps_ad)
  $SET(c_eps_S)
  $SET(kappa_ad)
  $SET(kappa_S)
  $SET(Omega_rot)
  $SET(T)

  !****

  function delta_p_(this) result(delta_p)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_p

    real(WP)              :: int_f
    integer               :: s
    real(WP), allocatable :: V_2(:)
    real(WP), allocatable :: c_1(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: f(:)

    ! Calculate the p-mode (large) frequency separation

    int_f = 0._WP

    seg_loop : do s = 1, this%n_s

       associate (ms => this%ms(s))

         V_2 = ms%V_2(ms%x)
         c_1 = ms%c_1(ms%x)
         Gamma_1 = ms%Gamma_1(ms%x)
         
         f = Gamma_1/(c_1*V_2)

         int_f = int_f + integrate(ms%x, f)

       end associate

    end do seg_loop

    delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/int_f
       
    ! Finish

    return

  end function delta_p_

  !****

  function delta_g_(this) result(delta_g)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_g

    real(WP)              :: int_f
    integer               :: s
    real(WP), allocatable :: As(:)
    real(WP), allocatable :: c_1(:)
    real(WP), allocatable :: f(:)

    ! Calculate the g-mode inverse period separation

    int_f = 0._WP

    seg_loop : do s = 1, this%n_s

       associate (ms => this%ms(s))

         As = ms%As(ms%x)
         c_1 = ms%c_1(ms%x)

         where(ms%x /= 0._WP)
            f = MAX(As/c_1, 0._WP))/ms%x
         elsewhere
            f = 0._WP
         end where

         int_f = int_f + integrate(ms%x, f)

       end associate

    end do seg_loop

    delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)*int_f/PI**2

    ! Finish

    return

  end function delta_g_

end module gyre_evol_model
