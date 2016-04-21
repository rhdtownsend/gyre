! Module   : gyre_evol_model
! Purpose  : stellar evolutionary model
!
! Copyright 2013-2016 Rich Townsend
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
  use gyre_grid
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_point
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $SET_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME}
  $endsub

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  type, extends (model_t) :: evol_model_t
     private
     type(grid_t)                  :: gr
     type(evol_seg_t), allocatable :: es(:)
     integer                       :: s_i
     integer                       :: s_o
     integer                       :: n_k
     real(WP), public              :: M_star
     real(WP), public              :: R_star
     real(WP), public              :: L_star
     logical                       :: add_center
     logical                       :: repair_As
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
     $SET_DECL(kap_ad)
     $SET_DECL(kap_S)
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
     $PROC_DECL(kap_ad)
     $PROC_DECL(kap_S)
     $PROC_DECL(Omega_rot)
     $PROC_DECL(dOmega_rot)
     $PROC_DECL(M_r)
     generic, public   :: M_r => M_r_1_, M_r_v_
     $PROC_DECL(P)
     generic, public   :: P => P_1_, P_v_
     $PROC_DECL(rho)
     generic, public   :: rho => rho_1_, rho_v_
     $PROC_DECL(T)
     generic, public   :: T => T_1_, T_v_
     procedure, public :: grid
     procedure, public :: vacuum
     procedure, public :: delta_p
     procedure, public :: delta_g
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

  function evol_model_t_ (x, M_star, R_star, L_star, ml_p) result (ml)

    real(WP), intent(in)          :: x(:)
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    type(model_par_t), intent(in) :: ml_p
    type(evol_model_t)            :: ml

    integer :: s

    ! Construct the evol_model_t

    ! Create the grid
    
    if (ml_p%add_center) then

       if (x(1) /= 0._WP) then

          ml%gr = grid_t([0._WP,x])
          ml%add_center = .TRUE.

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 100) 'Added central point'
100          format(3X,A)
          endif

       else

          ml%gr = grid_t(x)
          ml%add_center = .FALSE.

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 100) 'No need to add central point'
          endif

       endif

    else

       ml%gr = grid_t(x)
       ml%add_center = .FALSE.

    endif

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    ml%n_k = ml%gr%n_k

    ! Create segments

    allocate(ml%es(ml%s_i:ml%s_o))

    seg_loop : do s = ml%s_i, ml%s_o
       ml%es(s) = evol_seg_t(ml_p)
    end do seg_loop

    ml%repair_As = ml_p%repair_As

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ! Finish

    return

  end function evol_model_t_

  !****

  $define $SET $sub

  $local $NAME $1
  $local $F_C $2

  subroutine set_${NAME} (this, f)

    class(evol_model_t), intent(inout) :: this
    real(WP), intent(in)               :: f(:)

    real(WP), allocatable :: f_(:)
    integer               :: s
    integer               :: k_i
    integer               :: k_o

    ! Set the data for $NAME

    if (this%add_center) then
       f_ = [$F_C,f]
    else
       f_ = f
    endif

    $CHECK_BOUNDS(SIZE(f_),this%gr%n_k)

    seg_loop : do s = this%s_i, this%s_o

       k_i = this%gr%k_i(s)
       k_o = this%gr%k_o(s)

       associate (pt => this%gr%pt)
         call this%es(s)%set_${NAME}(pt(k_i:k_o)%x, f_(k_i:k_o))
       end associate
       
    end do seg_loop

    ! Finish

    return

  contains

    function f_c_ () result (f_c)

      real(WP) :: f_c

      ! Interpolate f at x=0 using parabolic fitting

      associate (x_1 => this%gr%pt(2)%x, &
                 x_2 => this%gr%pt(3)%x, &
                 f_1 => f(1), &
                 f_2 => f(2))

        f_c = (x_2**2*f_1 - x_1**2*f_2)/(x_2**2 - x_1**2)

      end associate

      ! Finish

      return

    end function f_c_

  end subroutine set_${NAME}

  $endsub

  $SET(V_2,f_c_())
  $SET(U,3._WP)
  $SET(c_1,f_c_())
  $SET(Gamma_1,f_c_())
  $SET(delta,f_c_())
  $SET(nabla_ad,f_c_())
  $SET(nabla,f_c_())
  $SET(beta_rad,f_c_())
  $SET(c_rad,f_c_())
  $SET(c_thm,f_c_())
  $SET(c_dif,f_c_())
  $SET(c_eps_ad,f_c_())
  $SET(c_eps_S,f_c_())
  $SET(kap_ad,f_c_())
  $SET(kap_S,f_c_())
  $SET(Omega_rot,f_c_())

  !****

  subroutine set_As (this, f)

    class(evol_model_t), intent(inout) :: this
    real(WP), intent(in)               :: f(:)

    real(WP), allocatable :: f_(:)
    integer               :: s
    integer               :: k_i
    integer               :: k_o

    ! Set the data for As

    if (this%add_center) then
       f_ = [0._WP,f]
    else
       f_ = f
    endif

    seg_loop : do s = this%s_i, this%s_o

       k_i = this%gr%k_i(s)
       k_o = this%gr%k_o(s)

       associate (pt => this%gr%pt)

         if (this%repair_As) then

            ! Repair the segment boundaries

            if (s > this%s_i .AND. k_i + 2 <= k_o) then
               f_(k_i) = f(k_i+1) + (pt(k_i)%x - pt(k_i+1)%x)*(f(k_i+2) - f(k_i+1))/(pt(k_i+2)%x - pt(k_i+1)%x)
            endif
               
            if (s < this%s_o .AND. k_o - 2 >= k_i) then
               f_(k_o) = f(k_o-1) + (pt(k_o)%x - pt(k_o-1)%x)*(f(k_o-1) - f(k_o-2))/(pt(k_o-1)%x - pt(k_o-2)%x)
            endif

         endif
               
         call this%es(s)%set_As(pt(k_i:k_o)%x, f_(k_i:k_o))

       end associate

    end do seg_loop

    ! Finish

    return

  end subroutine set_As

  !****

  $define $PROC_1 $sub

  $local $NAME $1

  function ${NAME}_1_ (this, pt) result (${NAME})

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: $NAME

    $ASSERT_DEBUG(pt%s >= this%s_i,Invalid segment)
    $ASSERT_DEBUG(pt%s <= this%s_o,Invalid segment)

    ! Evaluate $NAME at point pt

    $NAME = this%es(pt%s)%${NAME}(pt%x)

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
  $PROC_1(kap_ad)
  $PROC_1(kap_S)
  $PROC_1(Omega_rot)
  $PROC_1(dOmega_rot)

  !****

  function M_r_1_ (this, pt) result (M_r)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: M_r

    ! Evaluate the fractional mass coordinate at point pt

    M_r = this%M_star*(pt%x**3/this%c_1(pt))

    ! Finish

    return

  end function M_r_1_
    
  !****

  function P_1_ (this, pt) result (P)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: P

    ! Evaluate the total pressure at point pt

    P = (G_GRAVITY*this%M_star**2/(4._WP*PI*this%R_star**4))*&
        (this%U(pt)/(this%c_1(pt)**2*this%V_2(pt)))

    ! Finish

    return

  end function P_1_
    
  !****

  function rho_1_ (this, pt) result (rho)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: rho

    ! Evaluate the density at point pt

    rho = (this%M_star/(4._WP*PI*this%R_star**3))*(this%U(pt)/this%c_1(pt))

    ! Finish

    return

  end function rho_1_
    
  !****

  function T_1_ (this, pt) result (T)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: T

    ! Evaluate the temperature at point pt

    T = (3._WP*this%beta_rad(pt)*this%P(pt)/A_RADIATION)**0.25_WP

    ! Finish

    return

  end function T_1_
    
  !****

  $define $PROC_V $sub

  $local $NAME $1

  function ${NAME}_v_ (this, pt) result (${NAME})

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt(:)
    real(WP)                        :: ${NAME}(SIZE(pt))

    integer :: j

    ! Evaluate $NAME at points pt

    !$OMP PARALLEL DO
    do j = 1, SIZE(pt)
       ${NAME}(j) = this%${NAME}(pt(j))
    end do

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
  $PROC_V(kap_ad)
  $PROC_V(kap_S)
  $PROC_V(Omega_rot)
  $PROC_V(dOmega_rot)
  $PROC_V(M_r)
  $PROC_V(P)
  $PROC_V(rho)
  $PROC_V(T)

  !****

  function grid (this) result (gr)

    class(evol_model_t), intent(in) :: this
    type(grid_t)                    :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

  !****

  function vacuum (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    logical                         :: vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i,Invalid segment)
    $ASSERT_DEBUG(pt%s <= this%s_o,Invalid segment)

    ! Evaluate the vacuum condition at point pt

    vacuum = this%es(pt%s)%vacuum(pt%x)

    ! Finish

    return

  end function vacuum

  !****

  function delta_p (this)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_p

    real(WP) :: V_2(this%n_k)
    real(WP) :: c_1(this%n_k)
    real(WP) :: Gamma_1(this%n_k)
    real(WP) :: f(this%n_k)

    ! Calculate the p-mode (large) frequency separation

    associate (pt => this%gr%pt)

      V_2 = this%V_2(pt)
      c_1 = this%c_1(pt)
      Gamma_1 = this%Gamma_1(pt)

      f = Gamma_1/(c_1*V_2)

      $if ($GFORTRAN_PR_49636)
      delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/ &
                integrate(this%gr%pt%x, f)
      $else
      delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/ &
                integrate(pt%x, f)
      $endif

    end associate

    ! Finish

    return

  end function delta_p

  !****

  function delta_g (this)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_g

    real(WP) :: As(this%n_k)
    real(WP) :: c_1(this%n_k)
    real(WP) :: f(this%n_k)

    ! Calculate the g-mode inverse period separation

    associate (pt => this%gr%pt)

      As = this%As(pt)
      c_1 = this%c_1(pt)

      $if ($GFORTRAN_PR_49636)
      where (this%gr%pt%x /= 0._WP)
         f = SQRT(MAX(As/c_1, 0._WP))/this%gr%pt%x
      elsewhere
         f = 0._WP
      end where

      delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/PI**2* &
           integrate(this%gr%pt%x, f)
      $else
      where (pt%x /= 0._WP)
         f = SQRT(MAX(As/c_1, 0._WP))/pt%x
      elsewhere
         f = 0._WP
      end where

      delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/PI**2* &
                integrate(pt%x, f)
      $endif

    end associate

    ! Finish

    return

  end function delta_g

end module gyre_evol_model
