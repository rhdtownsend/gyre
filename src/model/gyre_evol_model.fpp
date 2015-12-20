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
  use gyre_evol_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $SET_DECL $sub
    $local $NAME $1
    procedure, public :: set_${NAME} => set_${NAME}_
  $endsub

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  type, extends (model_t) :: evol_model_t
     private
     integer, allocatable          :: s(:)
     real(WP), allocatable         :: x(:)
     type(evol_seg_t), allocatable :: es(:)
     real(WP), public              :: M_star
     real(WP), public              :: R_star
     real(WP), public              :: L_star
     integer                       :: n_k
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

  function evol_model_t_ (x, M_star, R_star, L_star, ml_p) result (ml)

    real(WP), intent(in)          :: x(:)
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    type(model_par_t), intent(in) :: ml_p
    type(evol_model_t)            :: ml

    integer :: s

    ! Construct the evol_model_t

    ml%n_k = SIZE(x)
    ml%n_s = ml%s(ml%n_k)

    ml%x = x
    ml%s = seg_indices_(x)

    ml%x_i = x(1)
    ml%x_o = x(ml%n_k)

    allocate(ml%es(ml%n_s))

    seg_loop : do s = 1, ml%n_s
       ml%es(s) = evol_seg_t(PACK(x, MASK=ml%s == s), ml_p)
    end do seg_loop

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ! Finish

    return

  contains

    function seg_indices_ (x) result (s)

      real(WP), intent(in) :: x(:)
      integer              :: s(SIZE(x))

      integer :: k

      $CHECK_BOUNDS(SIZE(mp),SIZE(x))

      ! Partition the array x into strictly-monotonic-increasing
      ! segments, by splitting at double points; return the resulting
      ! segment index of each element in s

      s(1) = 1

      x_loop : do k = 2, SIZE(x)
       
         if (x(k) == x(k-1)) then
            s(k) = s(k-1) + 1
         else
            s(k) = s(k-1)
         endif

      end do x_loop

      ! Finish

      return

    end function seg_indices_

  end function evol_model_t_

  !****

  $define $SET $sub

  $local $NAME $1

  subroutine set_${NAME}_ (this, f)

    class(evol_model_t), intent(inout) :: this
    real(WP), intent(in)               :: f(:)

    integer :: s

    $CHECK_BOUNDS(SIZE(f),this%n_k)

    ! Set the data for $NAME

    seg_loop : do s = 1, this%n_s
       call this%es(s)%set_${NAME}(PACK(f, MASK=this%s == s))
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

  function ${NAME}_1_ (this, s, x) result (${NAME})

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: $NAME

    ! Evaluate $NAME

    $NAME = this%es(s)%${NAME}(x)

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

  function M_r_1_ (this, s, x) result (M_r)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: M_r

    ! Evaluate the fractional mass coordinate

    M_r = this%M_star*(x**3*/this%c_1(s, x))

    ! Finish

    return

  end function M_r_1_
    
  !****

  function P_1_ (this, s, x) result (P)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: P

    ! Evaluate the total pressure

    P = (G_GRAVITY*this%M_star/(4._WP*PI*this%R_star**4))*&
        (x**2*this%U(s, x)/(this%c_1(s, x)*this%V(s, x)))

    ! Finish

    return

  end function P_1_
    
  !****

  function rho_1_ (this, s, x) result (rho)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: rho

    ! Evaluate the density

    rho = (this%M_star/(4._WP*PI*this%R_star)**3*(this%U(s, x)/this%c_1(s, x))

    ! Finish

    return

  end function rho_1_
    
  !****

  function T_1_ (this, s, x) result (T)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: T

    ! Evaluate the temperature

    T = (3._WP*this%beta_rad(s, x)*this%P(s, x)/A_RADIATION)**0.25_WP

    ! Finish

    return

  end function T_1_
    
  !****

  $define $PROC_V $sub

  $local $NAME $1

  function ${NAME}_v_ (this, co) result (${NAME})

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: s(:)
    real(WP), intent(in)            :: x(:)
    real(WP)                        :: ${NAME}(SIZE(s))

    integer :: k

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Evaluate $NAME

    !$OMP PARALLEL DO
    do k = 1, SIZE(s)
       ${NAME}(k) = this%${NAME}(s(k), x(k))
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
  $PROC_V(kappa_ad)
  $PROC_V(kappa_S)
  $PROC_V(Omega_rot)
  $PROC_V(dOmega_rot)
  $PROC_V(M_r)
  $PROC_V(P)
  $PROC_V(rho)
  $PROC_V(T)

  !****

  subroutine scaffold_ (this, s, x)

    class(evol_model_t), intent(in)    :: this
    integer, allocatable, intent(out)  :: s(:)
    real(WP), allocatable, intent(out) :: x(:)

    ! Return the grid scaffold

    s = this%s
    x = this%x

    ! Finish

    return

  end function scaffold_

  !****

  function delta_p_ (this) result (delta_p)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_p

    real(WP) :: V_2(this%n_k)
    real(WP) :: c_1(this%n_k)
    real(WP) :: Gamma_1(this%n_k)
    real(WP) :: f(this%n_k)

    ! Calculate the p-mode (large) frequency separation

    associate (s => this%s, &
               x => this%x)

      V_2 = this%V_2(s, x)
      c_1 = this%c_1(s, x)

      Gamma_1 = this%Gamma_1(s, x)

      f = Gamma_1/(c_1*V_2)

      delta_p = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/ &
                integrate(x, f)

    end associate
       
    ! Finish

    return

  end function delta_p_

  !****

  function delta_g_ (this) result (delta_g)

    class(evol_model_t), intent(in) :: this
    real(WP)                        :: delta_g

    real(WP) :: As(this%n_k)
    real(WP) :: c_1(this%n_k)
    real(WP) :: f(this%n_k)

    ! Calculate the g-mode inverse period separation

    associate (s => this%s, &
               x => this%x)

      As = this%As(s, x)
      c_1 = this%c_1(s, x)

      where (x /= 0._WP)
         f = MAX(As/c_1, 0._WP))/x
      elsewhere
         f = 0._WP
      end where

      delta_g = 0.5_WP*SQRT(G_GRAVITY*this%M_star/this%R_star**3)/PI**2* &
           integrate(x, f)

    end associate

    ! Finish

    return

  end function delta_g_

end module gyre_evol_model
