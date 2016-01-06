! Module   : gyre_poly_model
! Purpose  : stellar polytropic model
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

module gyre_poly_model

  ! Uses

  use core_kinds

  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_poly_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  type, extends (model_t) :: poly_model_t
     private
     type(poly_seg_t), allocatable :: ps(:)
     real(WP), allocatable         :: x(:)
     integer, allocatable          :: k_i(:)
     integer, allocatable          :: k_o(:)
     integer                       :: n_k
   contains
     private
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
     procedure, public :: x_i
     procedure, public :: x_o
     procedure, public :: x_s
  end type poly_model_t

  ! Interfaces

  interface poly_model_t
     module procedure poly_model_t_
  end interface poly_model_t

  ! Access specifiers

  private

  public :: poly_model_t

  ! Procedures

contains

  function poly_model_t_ (xi, Theta, dTheta, n_poly, Gamma_1, ml_p) result (ml)

    real(WP), intent(in)          :: xi(:)
    real(WP), intent(in)          :: Theta(:)
    real(WP), intent(in)          :: dTheta(:)
    real(WP), intent(in)          :: n_poly(:)
    real(WP), intent(in)          :: Gamma_1
    type(model_par_t), intent(in) :: ml_p
    type(poly_model_t)            :: ml

    real(WP), allocatable :: mu(:)
    real(WP), allocatable :: B(:)
    integer               :: s

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    ! Construct the poly_model_t

    ml%n_k = SIZE(xi)

    ml%x = xi/xi(ml%n_k)
    call seg_indices(ml%x, ml%k_i, ml%k_o)

    ml%n_s = SIZE(ml%k_i)

    $CHECK_BOUNDS(SIZE(n_poly),ml%n_s)

    allocate(mu(ml%n_s+1))
    allocate(B(ml%n_s))

    call eval_mu_B_(xi, dTheta, ml%k_i, ml%k_o, mu, B)

    allocate(ml%ps(ml%n_s))

    seg_loop : do s = 1, ml%n_s

       associate (k_i => ml%k_i(s), &
                  k_o => ml%k_o(s))

         ml%ps(s) = poly_seg_t(xi(k_i:k_o), Theta(k_i:k_o), dTheta(k_i:k_o), &
                               mu(s), mu(ml%n_s+1), B(s), xi(ml%n_k), n_poly(s), Gamma_1)

       end associate

    end do seg_loop

    ! Finish

    return

  end function poly_model_t_

  !****

  subroutine eval_mu_B_ (xi, dTheta, k_i, k_o, mu, B)

    real(WP), intent(in)  :: xi(:)
    real(WP), intent(in)  :: dTheta(:)
    integer, intent(in)   :: k_i(:)
    integer, intent(in)   :: k_o(:)
    real(WP), intent(out) :: mu(:)
    real(WP), intent(out) :: B(:)

    integer  :: n_s
    integer  :: s
    real(WP) :: v_i
    real(WP) :: v_o

    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    $CHECK_BOUNDS(SIZE(k_o),SIZE(k_i))

    $CHECK_BOUNDS(SIZE(mu),SIZE(k_i)+1)
    $CHECK_BOUNDS(SIZE(B),SIZE(k_i))

    ! Calculate the mass coordinate mu and jump factor B at segment
    ! boundaries

    n_s = SIZE(k_i)

    mu(1) = 0._WP
    B(1) = 1._WP
    
    seg_loop : do s = 1, n_s-1

       v_i = xi(k_i(s))**2*dTheta(k_i(s))
       v_o = xi(k_o(s))**2*dTheta(k_o(s))

       mu(s+1) = mu(s) - (v_o - v_i)/B(s)
       B(s+1) = B(s)*dTheta(k_i(s+1))/dTheta(k_o(s))

    end do seg_loop

    v_i = xi(k_i(s))**2*dTheta(k_i(s))
    v_o = xi(k_o(s))**2*dTheta(k_o(s))

    mu(s+1) = mu(s) - (v_o - v_i)/B(s)
    
    ! Finish

    return

  end subroutine eval_mu_B_

  !****

  $define $PROC_1 $sub

  $local $NAME $1

  function ${NAME}_1_ (this, s, x) result (${NAME})

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: $NAME

    $ASSERT_DEBUG(s >= 1,Invalid segment index)
    $ASSERT_DEBUG(s <= this%n_s,Invalid segment index)

    ! Evaluate $NAME

    $NAME = this%ps(s)%${NAME}(x)

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
  $PROC_1(Omega_rot)
  $PROC_1(dOmega_rot)

  !****

  $define $PROC_1_NULL $sub

  $local $NAME $1

  function ${NAME}_1_ (this, s, x) result (${NAME})

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), intent(in)            :: x
    real(WP)                        :: $NAME

    $ABORT(Polytropic model does not define $NAME)

    ! (This line to prevent unset warnings)

    $NAME = 0._WP

    ! Finish

    return

  end function ${NAME}_1_

  $endsub

  $PROC_1_NULL(nabla)
  $PROC_1_NULL(beta_rad)
  $PROC_1_NULL(c_rad)
  $PROC_1_NULL(dc_rad)
  $PROC_1_NULL(c_thm)
  $PROC_1_NULL(c_dif)
  $PROC_1_NULL(c_eps_ad)
  $PROC_1_NULL(c_eps_S)
  $PROC_1_NULL(kappa_ad)
  $PROC_1_NULL(kappa_S)

  !****

  $define $PROC_V $sub

  $local $NAME $1

  function ${NAME}_v_ (this, s, x) result (${NAME})

    class(poly_model_t), intent(in) :: this
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

  !****

  function x_i (this, s)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP)                        :: x_i

    $ASSERT_DEBUG(s >= 1,Invalid segment index)
    $ASSERT_DEBUG(s <= this%n_s,Invalid segment index)

    ! Return the inner x value for segment s

    x_i = this%x(this%k_i(s))

    ! Finish

    return

  end function x_i

  !****

  function x_o (this, s)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP)                        :: x_o

    $ASSERT_DEBUG(s >= 1,Invalid segment index)
    $ASSERT_DEBUG(s <= this%n_s,Invalid segment index)

    ! Return the outer x value for segment s

    x_o = this%x(this%k_o(s))

    ! Finish

    return

  end function x_o

  !****

  function x_s (this, s)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: s
    real(WP), allocatable           :: x_s(:)

    $ASSERT_DEBUG(s >= 1,Invalid segment index)
    $ASSERT_DEBUG(s <= this%n_s,Invalid segment index)

    ! Return the model grid for segment s

    x_s = this%x(this%k_i(s):this%k_o(s))

    ! Finish

    return

  end function x_s

end module gyre_poly_model
