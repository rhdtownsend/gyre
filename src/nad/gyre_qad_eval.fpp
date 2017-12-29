! Module   : gyre_qad_eval
! Purpose  : quasiadiabatic eigenfunction evaluation
!
! Copyright 2017 Rich Townsend
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

module gyre_qad_eval

  ! Uses

  use core_kinds

  use gyre_calc
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_nad_eqns
  use gyre_osc_par
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: qad_eval_t
     private
     type(context_t), pointer :: cx => null()
     type(nad_eqns_t)         :: eq
     type(grid_t)             :: gr
     integer, public          :: n_k
   contains
     private
     procedure, public :: y_qad
  end type qad_eval_t

  ! Interfaces

  interface qad_eval_t
     module procedure qad_eval_t_
  end interface qad_eval_t

  ! Access specifiers

  private

  public :: qad_eval_t

  ! Procedures

contains

  function qad_eval_t_ (cx, gr, md_p, os_p) result (qe)

    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(qad_eval_t)                     :: qe

    ! Construct the qad_eval_t

    ! Initialize the equations

    qe%eq = nad_eqns_t(cx, md_p, os_p)
    call qe%eq%stencil(gr%pt)

    ! Other initializations

    qe%cx => cx
    qe%gr = gr

    qe%n_k = gr%n_k

    ! Finish

    return

  end function qad_eval_t_

  !****

  function y_qad (this, st, y_ad)

    class(qad_eval_t), intent(inout) :: this
    class(c_state_t), intent(in)     :: st
    real(WP), intent(in)             :: y_ad(:,:)
    complex(WP)                      :: y_qad(6,this%n_k)

    integer         :: k
    integer         :: s
    complex(WP)     :: xA(6,6)
    complex(WP)     :: xA_5(6,this%n_k)
    complex(WP)     :: xA_6(6,this%n_k)
    complex(WP)     :: dy_6(this%n_k)

    $CHECK_BOUNDS(SIZE(y_ad, 1),4)
    $CHECK_BOUNDS(SIZE(y_ad, 2),this%n_k)

    ! Construct quasi-adiabatic eigenfunctions y_qad from adiabatic
    ! eigenfrequency omega_ad and eigenfunctions y_ad

    ! Copy over adiabatic eigenfunctions

    y_qad(1:4,:) = y_ad

    ! Evaluate components of the non-adiabatic RHS matrix
    ! corresponding to the energy conservation and transport equations

    !$OMP PARALLEL DO PRIVATE (xA)
    do k = 1, this%n_k

       xA = this%eq%xA(k, st)

       xA_5(:,k) = xA(5,:)
       xA_6(:,k) = xA(6,:)

    end do

    ! Evaluate the luminosity perturbation eigenfunction

    where (this%gr%pt%x /= 0._WP)
       y_qad(6,:) = -(xA_5(1,:)*y_qad(1,:) + xA_5(2,:)*y_qad(2,:) + xA_5(3,:)*y_qad(3,:) + xA_5(4,:)*y_qad(4,:))/xA_5(6,:)
    elsewhere
       y_qad(6,:) = 0._WP
    end where

    ! Evaluate the gradient of the luminosity perturbation
    ! eigenfunction, segment-by-segment
    
    seg_loop : do s = this%gr%s_i(), this%gr%s_o()
       associate (k_i => this%gr%k_s_i(s), k_o => this%gr%k_s_o(s))
         dy_6(k_i:k_o) = this%gr%pt(k_i:k_o)%x*deriv(this%gr%pt(k_i:k_o)%x, REAL(y_qad(6,k_i:k_o), WP), 'MONO')
       end associate
    end do seg_loop

    ! Evaluate the entropy perturbation eigenfunction

    where (this%gr%pt%x /= 0._WP)
       y_qad(5,:) = (dy_6 - (xA_6(1,:)*y_qad(1,:) + xA_6(2,:)*y_qad(2,:) + xA_6(3,:)*y_qad(3,:) + &
                             xA_6(4,:)*y_qad(4,:) + xA_6(6,:)*y_qad(6,:)))/xA_6(5,:)
    elsewhere
       y_qad(5,:) = 0._WP
    end where

    ! Finish

    return

  end function y_qad

end module gyre_qad_eval
