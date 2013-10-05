! Module   : gyre_ivp_magnus
! Purpose  : solve initial-value problems (interface, for Magnus schemes)
!
! Copyright 2013 Rich Townsend
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

module gyre_ivp_magnus

  ! Uses

  use core_kinds
  use core_constants

  use gyre_jacobian
  use gyre_ext_arith
  use gyre_ivp
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract, extends (ivp_t) :: ivp_magnus_t
   contains
     private
     procedure, public                          :: solve
     procedure, public                          :: recon
     procedure(eval_dOmega_i), deferred, public :: eval_dOmega
  end type ivp_magnus_t

  ! Interfaces

  abstract interface

     subroutine eval_dOmega_i (this, omega, x_a, x_b, dOmega)
       use core_kinds
       import ivp_magnus_t
       class(ivp_magnus_t), intent(in) :: this
       complex(WP), intent(in)         :: omega
       real(WP), intent(in)            :: x_a
       real(WP), intent(in)            :: x_b
       complex(WP), intent(out)        :: dOmega(:,:)
     end subroutine eval_dOmega_i

  end interface

  ! Access specifiers

  private

  public :: ivp_magnus_t

  ! Procedures

contains

  subroutine solve (this, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(ivp_magnus_t), intent(in)  :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    logical :: UPWIND = .TRUE.

    logical     :: use_real_
    complex(WP) :: dOmega(this%n_e,this%n_e)
    real(WP)    :: dOmega_r(this%n_e,this%n_e)
    complex(WP) :: lambda(this%n_e)
    complex(WP) :: V_l(this%n_e,this%n_e)
    complex(WP) :: V_r(this%n_e,this%n_e)
    real(WP)    :: dx
    integer     :: i
    complex(WP) :: V_pos(this%n_e,this%n_e)
    complex(WP) :: V_neg(this%n_e,this%n_e)

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    ! Evaluate the Magnus slope matrix

    call this%eval_dOmega(omega, x_a, x_b, dOmega)

    ! Decompose it

    if(use_real_) then
       dOmega_r = REAL(dOmega)
       call eigen_decompose(dOmega_r, lambda, V_l, V_r)
    else
       call eigen_decompose(dOmega, lambda, V_l, V_r)
    endif

    ! Set up the solution matrices and scales, using 'upwinding' for stability

    dx = x_b - x_a

    $block

    $if($DOUBLE_PRECISION)
    $local $X Z
    $else
    $local $X C
    $endif

    if(UPWIND) then

       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_pos(1,i), 1)
          if(REAL(lambda(i)) >= 0._WP) then
             call ${X}SCAL(this%n_e, EXP(-lambda(i)*dx), V_pos(1,i), 1)
          endif
       end do
    
       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          if(REAL(lambda(i)) < 0._WP) then
             call ${X}SCAL(this%n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
          endif
       end do

       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(-1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l, this%n_e)

       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_pos, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_r, this%n_e)

       S = exp(ext_complex(SUM(lambda, MASK=REAL(lambda) >= 0._WP)*dx))

    else

       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          call ${X}SCAL(this%n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
       end do
    
       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l, this%n_e)

       do i = 1,this%n_e
          E_r(:,i) = 0._WP
          E_r(i,i) = -1._WP
       end do

       S = ext_complex(1._WP)

    endif
    
    $endblock

    ! Finish

    return

  end subroutine solve

!****

  subroutine recon (this, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(ivp_magnus_t), intent(in) :: this
    complex(WP), intent(in)         :: omega
    real(WP), intent(in)            :: x_a
    real(WP), intent(in)            :: x_b
    complex(WP), intent(in)         :: y_a(:)
    complex(WP), intent(in)         :: y_b(:)
    real(WP), intent(in)            :: x(:)
    complex(WP), intent(out)        :: y(:,:)
    logical, intent(in), optional   :: use_real

    logical     :: use_real_
    complex(WP) :: dOmega(this%n_e,this%n_e)
    real(WP)    :: dOmega_r(this%n_e,this%n_e)
    complex(WP) :: lambda(this%n_e)
    complex(WP) :: V_l(this%n_e,this%n_e)
    complex(WP) :: V_r(this%n_e,this%n_e)
    integer     :: i
    complex(WP) :: exp_a(this%n_e)
    complex(WP) :: exp_b(this%n_e)

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    ! Evaluate the Magnus slope matrix

    call this%eval_dOmega(omega, x_a, x_b, dOmega)

    ! Decompose it

    if(use_real_) then
       dOmega_r = REAL(dOmega)
       call eigen_decompose(dOmega_r, lambda, V_l, V_r)
    else
       call eigen_decompose(dOmega, lambda, V_l, V_r)
    endif

    ! Do the stabilized (both-boundaries) Magnus reconstruction

    recon_loop : do i = 1,SIZE(x)

       exp_a = MERGE(EXP(lambda*(x(i) - x_a)), CMPLX(0._WP, KIND=WP), REAL(lambda) < 0._WP .EQV. x_b > x_a)
       exp_b = MERGE(EXP(lambda*(x(i) - x_b)), CMPLX(0._WP, KIND=WP), REAL(lambda) >= 0._WP .EQV. x_b > x_a)
       
       y(:,i) = MATMUL(V_r, MATMUL(diagonal_matrix(exp_a), MATMUL(V_l, y_a)) + &
                            MATMUL(diagonal_matrix(exp_b), MATMUL(V_l, y_b)))

    end do recon_loop

    ! Finish

    return

  end subroutine recon

end module gyre_ivp_magnus
