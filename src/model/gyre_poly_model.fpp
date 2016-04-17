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

  use gyre_grid
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_point
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
     type(grid_t)                  :: gr
     type(poly_seg_t), allocatable :: ps(:)
     integer                       :: s_i
     integer                       :: s_o
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
     $PROC_DECL(kap_ad)
     $PROC_DECL(kap_S)
     $PROC_DECL(Omega_rot)
     $PROC_DECL(dOmega_rot)
     procedure, public :: grid
     procedure, public :: vacuum
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

  function poly_model_t_ (xi, Theta, dTheta, n_poly, Delta_d, Gamma_1, ml_p) result (ml)

    real(WP), intent(in)          :: xi(:)
    real(WP), intent(in)          :: Theta(:)
    real(WP), intent(in)          :: dTheta(:)
    real(WP), intent(in)          :: n_poly(:)
    real(WP), intent(in)          :: Delta_d(:)
    real(WP), intent(in)          :: Gamma_1
    type(model_par_t), intent(in) :: ml_p
    type(poly_model_t)            :: ml

    integer  :: n_k
    real(WP) :: x(SIZE(xi))
    integer  :: k_i(SIZE(n_poly))
    integer  :: k_o(SIZE(n_poly))
    real(WP) :: mu(SIZE(n_poly)+1)
    real(WP) :: B(SIZE(n_poly))
    real(WP) :: t(SIZE(n_poly))
    integer  :: s
    integer  :: i

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    $CHECK_BOUNDS(SIZE(Delta_d),SIZE(n_poly)-1)

    ! Construct the poly_model_t

    n_k = SIZE(xi)

    x = xi/xi(n_k)

    ! Create the grid

    ml%gr = grid_t(x)

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    ml%n_k = ml%gr%n_k

    $CHECK_BOUNDS(SIZE(n_poly),ml%s_o-ml%s_i+1)

    ! Evaluate mu, B and t

    do s = ml%s_i, ml%s_o
       i = s - ml%s_i + 1
       k_i(i) = ml%gr%k_i(s)
       k_o(i) = ml%gr%k_o(s)
    end do

    call eval_mu_B_t_(xi, Theta, dTheta, n_poly, Delta_d, k_i, k_o, mu, B, t)

    ! Create segments

    allocate(ml%ps(ml%s_i:ml%s_o))

    seg_loop : do s = ml%s_i, ml%s_o

       i = s - ml%s_i + 1

       ml%ps(s) = poly_seg_t(x(k_i(i):k_o(i)), Theta(k_i(i):k_o(i)), dTheta(k_i(i):k_o(i)), &
                             mu(i), mu(ml%s_o-ml%s_i+1), xi(n_k), &
                             n_poly(i), B(i), t(i), Gamma_1)

    end do seg_loop

    ! Finish

    return

  end function poly_model_t_

  !****

  subroutine eval_mu_B_t_(xi, Theta, dTheta, n_poly, Delta_d, k_i, k_o, mu, B, t)

    real(WP), intent(in)  :: xi(:)
    real(WP), intent(in)  :: Theta(:)
    real(WP), intent(in)  :: dTheta(:)
    real(WP), intent(in)  :: n_poly(:)
    real(WP), intent(in)  :: Delta_d(:)
    integer, intent(in)   :: k_i(:)
    integer, intent(in)   :: k_o(:)
    real(WP), intent(out) :: mu(:)
    real(WP), intent(out) :: B(:)
    real(WP), intent(out) :: t(:)

    integer  :: n_s
    integer  :: s
    real(WP) :: v_i
    real(WP) :: v_o

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    $CHECK_BOUNDS(SIZE(Delta_d),SIZE(n_poly)-1)

    $CHECK_BOUNDS(SIZE(k_i),SIZE(n_poly))
    $CHECK_BOUNDS(SIZE(k_o),SIZE(n_poly))

    $CHECK_BOUNDS(SIZE(mu),SIZE(n_poly)+1)
    $CHECK_BOUNDS(SIZE(B),SIZE(n_poly))
    $CHECK_BOUNDS(SIZE(t),SIZE(n_poly))

    ! Calculate the mass coordinate mu at segment boundaries, together
    ! with the B and t parameter (see mixed-poly-mod.pdf)

    n_s = SIZE(k_i)

    mu(1) = 0._WP
    B(1) = 1._WP
    t(1) = 1._WP

    seg_loop : do s = 1, n_s-1

       v_i = xi(k_i(s))**2*dTheta(k_i(s))
       v_o = xi(k_o(s))**2*dTheta(k_o(s))

       mu(s+1) = mu(s) - (v_o - v_i)*t(s)/B(s)

       t(s+1) = t(s)*EXP(n_poly(s)*LOG(Theta(k_o(s))) + Delta_d(s))

       B(s+1) = dTheta(k_i(s+1))/dTheta(k_o(s))*(t(s+1)/t(s))*B(s)

    end do seg_loop

    v_i = xi(k_i(s))**2*dTheta(k_i(s))
    v_o = xi(k_o(s))**2*dTheta(k_o(s))

    mu(s+1) = mu(s) - (v_o - v_i)*t(s)/B(s)

    ! Finish

    return

  end subroutine eval_mu_B_t_

  !****

  $define $PROC_1 $sub

  $local $NAME $1

  function ${NAME}_1_ (this, pt) result (${NAME})

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: $NAME

    $ASSERT_DEBUG(pt%s >= this%s_i,Invalid index)
    $ASSERT_DEBUG(pt%s <= this%s_o,Invalid index)

    ! Evaluate $NAME

    $NAME = this%ps(pt%s)%${NAME}(pt%x)

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

  function ${NAME}_1_ (this, pt) result (${NAME})

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
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
  $PROC_1_NULL(kap_ad)
  $PROC_1_NULL(kap_S)

  !****

  $define $PROC_V $sub

  $local $NAME $1

  function ${NAME}_v_ (this, pt) result (${NAME})

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt(:)
    real(WP)                        :: ${NAME}(SIZE(pt))

    integer :: j

    ! Evaluate $NAME

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

  !****

  function grid (this) result (gr)

    class(poly_model_t), intent(in) :: this
    type(grid_t)                    :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

  !****

  function vacuum (this, pt)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    logical                         :: vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i,Invalid segment)
    $ASSERT_DEBUG(pt%s <= this%s_o,Invalid segment)

    ! Evaluate the vacuum condition

    vacuum = this%ps(pt%s)%vacuum(pt%x)

    ! Finish

    return

  end function vacuum

end module gyre_poly_model
