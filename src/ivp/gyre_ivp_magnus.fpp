! Module   : gyre_ivp_magnus
! Purpose  : solve IVPs across fixed intervals using Magnus integrators
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
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: solve_magnus_GL2
  public :: solve_magnus_GL4
  public :: solve_magnus_GL6
  public :: recon_magnus_GL2
  public :: recon_magnus_GL4
  public :: recon_magnus_GL6
  public :: abscissa_magnus_GL2
  public :: abscissa_magnus_GL4
  public :: abscissa_magnus_GL6

  ! Procedures

contains

  subroutine solve_magnus_GL2 (jc, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using a 2nd-order
    ! Gauss-Legendre Magnus scheme

    call eval_dOmega_GL2(jc, omega, x_a, x_b, dOmega)
    call solve_magnus(dOmega, x_b-x_a, E_l, E_r, S, use_real)

    ! Finish
    
  end subroutine solve_magnus_GL2

!****
  
  subroutine solve_magnus_GL4 (jc, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using a 4th-order
    ! Gauss-Legendre Magnus scheme

    call eval_dOmega_GL4(jc, omega, x_a, x_b, dOmega)
    call solve_magnus(dOmega, x_b-x_a, E_l, E_r, S, use_real)

    ! Finish
    
  end subroutine solve_magnus_GL4
  
!****
  
  subroutine solve_magnus_GL6 (jc, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using a 6th-order
    ! Gauss-Legendre Magnus scheme

    call eval_dOmega_GL6(jc, omega, x_a, x_b, dOmega)
    call solve_magnus(dOmega, x_b-x_a, E_l, E_r, S, use_real)

    ! Finish
    
  end subroutine solve_magnus_GL6
  
!****

  subroutine solve_magnus (dOmega, dx, E_l, E_r, S, use_real)

    complex(WP), intent(inout)       :: dOmega(:,:)
    real(WP), intent(in)             :: dx
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    logical :: UPWIND = .TRUE.

    logical     :: use_real_
    real(WP)    :: dOmega_r(SIZE(dOmega, 1),SIZE(dOmega, 2))
    complex(WP) :: lambda(SIZE(dOmega, 1))
    complex(WP) :: V_l(SIZE(dOmega, 1),SIZE(dOmega, 1))
    complex(WP) :: V_r(SIZE(dOmega, 1),SIZE(dOmega, 1))
    integer     :: n_e
    integer     :: i
    complex(WP) :: V_pos(SIZE(dOmega, 1),SIZE(dOmega, 1))
    complex(WP) :: V_neg(SIZE(dOmega, 1),SIZE(dOmega, 1))

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    ! Decompose the Magnus slope matrix

    if(use_real_) then
       dOmega_r = REAL(dOmega)
       call eigen_decompose(dOmega_r, lambda, V_l, V_r)
    else
       call eigen_decompose(dOmega, lambda, V_l, V_r)
    endif

    ! Set up the solution matrices and scales, using 'upwinding' for stability

    n_e = SIZE(dOmega, 1)

    $block

    $if($DOUBLE_PRECISION)
    $local $X Z
    $else
    $local $X C
    $endif

    if(UPWIND) then

       do i = 1,n_e
          call ${X}COPY(n_e, V_r(1,i), 1, V_pos(1,i), 1)
          if(REAL(lambda(i)) >= 0._WP) then
             call ${X}SCAL(n_e, EXP(-lambda(i)*dx), V_pos(1,i), 1)
          endif
       end do
    
       do i = 1,n_e
          call ${X}COPY(n_e, V_r(1,i), 1, V_neg(1,i), 1)
          if(REAL(lambda(i)) < 0._WP) then
             call ${X}SCAL(n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
          endif
       end do

       call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(-1._WP, KIND=WP), &
                     V_neg, n_e, V_l, n_e, CMPLX(0._WP, KIND=WP), &
                     E_l, n_e)

       call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(1._WP, KIND=WP), &
                     V_pos, n_e, V_l, n_e, CMPLX(0._WP, KIND=WP), &
                     E_r, n_e)

       S = exp(ext_complex(SUM(lambda, MASK=REAL(lambda) >= 0._WP)*dx))

    else

       do i = 1,n_e
          call ${X}COPY(n_e, V_r(1,i), 1, V_neg(1,i), 1)
          call ${X}SCAL(n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
       end do
    
       call ${X}GEMM('N', 'N', n_e, n_e, n_e, CMPLX(1._WP, KIND=WP), &
                     V_neg, n_e, V_l, n_e, CMPLX(0._WP, KIND=WP), &
                     E_l, n_e)

       do i = 1,n_e
          E_r(:,i) = 0._WP
          E_r(i,i) = -1._WP
       end do

       S = ext_complex(1._WP)

    endif
    
    $endblock

    ! Finish

    return

  end subroutine solve_magnus

!****

  subroutine recon_magnus_GL2 (jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)
    logical, intent(in), optional :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(y_a),jc%n_e)
    $CHECK_BOUNDS(SIZE(y_b),jc%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the solution using a 2nd-order Gauss-Legendre Magnus
    ! scheme

    call eval_dOmega_GL2(jc, omega, x_a, x_b, dOmega)
    call recon_magnus (dOmega, x_a, x_b, y_a, y_b, x, y, use_real)

    ! Finish

    return

  end subroutine recon_magnus_GL2

!****

  subroutine recon_magnus_GL4 (jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)
    logical, intent(in), optional :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(y_a),jc%n_e)
    $CHECK_BOUNDS(SIZE(y_b),jc%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the solution using a 4th-order Gauss-Legendre Magnus
    ! scheme

    call eval_dOmega_GL4(jc, omega, x_a, x_b, dOmega)
    call recon_magnus (dOmega, x_a, x_b, y_a, y_b, x, y, use_real)

    ! Finish

    return

  end subroutine recon_magnus_GL4

!****

  subroutine recon_magnus_GL6 (jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)
    logical, intent(in), optional :: use_real

    complex(WP) :: dOmega(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(y_a),jc%n_e)
    $CHECK_BOUNDS(SIZE(y_b),jc%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the solution using a 6th-order Gauss-Legendre Magnus
    ! scheme

    call eval_dOmega_GL6(jc, omega, x_a, x_b, dOmega)
    call recon_magnus (dOmega, x_a, x_b, y_a, y_b, x, y, use_real)

    ! Finish

    return

  end subroutine recon_magnus_GL6

!****

  subroutine recon_magnus (dOmega, x_a, x_b, y_a, y_b, x, y, use_real)

    complex(WP), intent(inout)    :: dOmega(:,:)
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)
    logical, intent(in), optional :: use_real

    logical     :: use_real_
    real(WP)    :: dOmega_r(SIZE(dOmega, 1),SIZE(dOmega, 2))
    complex(WP) :: lambda(SIZE(dOmega, 1))
    complex(WP) :: V_l(SIZE(dOmega, 1),SIZE(dOmega, 1))
    complex(WP) :: V_r(SIZE(dOmega, 1),SIZE(dOmega, 1))
    integer     :: i
    complex(WP) :: exp_a(SIZE(dOmega,1))
    complex(WP) :: exp_b(SIZE(dOmega,1))

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    ! Reconstruct the solution within the interval x_a -> x_b using
    ! the Magnus slope matrix

    ! Decompose the matrix

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

  end subroutine recon_magnus

!****

  subroutine eval_dOmega_GL2 (jc, omega, x_a, x_b, dOmega)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(out)      :: dOmega(:,:)

    real(WP)    :: dx
    real(WP)    :: x(1)
    complex(WP) :: A(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(dOmega, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(dOmega, 2),jc%n_e)

    ! Evaluate the Magnus slope matrix for a 2nd-order Gauss-Legendre
    ! Magnus scheme

    ! Evaluate the Jacobian

    dx = x_b - x_a
    x = abscissa_magnus_GL2(x_a, x_b)

    call jc%eval(x(1), omega, A)

    ! Evaluate the slope matrix

    dOmega = A

    ! Finish

    return

  end subroutine eval_dOmega_GL2

!****

  subroutine eval_dOmega_GL4 (jc, omega, x_a, x_b, dOmega)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(out)      :: dOmega(:,:)

    real(WP)    :: dx
    real(WP)    :: x(2)
    complex(WP) :: A(jc%n_e,jc%n_e,2)
    complex(WP) :: dalpha(jc%n_e,jc%n_e,2)

    $CHECK_BOUNDS(SIZE(dOmega, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(dOmega, 2),jc%n_e)

    ! Evaluate the Magnus slope matrix for a 4th-order Gauss-Legendre
    ! Magnus scheme

    ! Evaluate the Jacobian

    dx = x_b - x_a
    x = abscissa_magnus_GL4(x_a, x_b)

    call jc%eval(x(1), omega, A(:,:,1))
    call jc%eval(x(2), omega, A(:,:,2))

    ! Evaluate the Magnus slope matrix (Blanes et al. 2009, eqns. 243
    ! and 253; note that the 12 in the denominator of their expression
    ! for alpha_2 is erroneous)

    dalpha(:,:,1) = 0.5_WP*(A(:,:,1) + A(:,:,2))
    dalpha(:,:,2) = SQRT(3._WP)*(A(:,:,2) - A(:,:,1))

    dOmega = dalpha(:,:,1) - dx*commutator(dalpha(:,:,1), dalpha(:,:,2))/12._WP
    
    ! Finish

    return

  end subroutine eval_dOmega_GL4

!****

  subroutine eval_dOmega_GL6 (jc, omega, x_a, x_b, dOmega)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(out)      :: dOmega(:,:)

    real(WP)    :: dx
    real(WP)    :: x(3)
    complex(WP) :: A(jc%n_e,jc%n_e,3)
    complex(WP) :: dalpha(jc%n_e,jc%n_e,3)
    complex(WP) :: dC(jc%n_e,jc%n_e,2)

    $CHECK_BOUNDS(SIZE(dOmega, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(dOmega, 2),jc%n_e)

    ! Evaluate the Magnus slope matrix for a 6th-order Gauss-Legendre
    ! Magnus scheme

    ! Evaluate Jacobians

    dx = x_b - x_a
    x = abscissa_magnus_GL6(x_a, x_b)

    call jc%eval(x(1), omega, A(:,:,1))
    call jc%eval(x(2), omega, A(:,:,2))
    call jc%eval(x(3), omega, A(:,:,3))

    ! Evaluate the Magnus slope matrix (Blanes et al. 2009, eqns. 251 and 257)

    dalpha(:,:,1) = A(:,:,2)
    dalpha(:,:,2) = SQRT(15._WP)*(A(:,:,3) - A(:,:,1))/3
    dalpha(:,:,3) = 10*(A(:,:,3) - 2*A(:,:,2) + A(:,:,1))/3

    dC(:,:,1) = dx*commutator(dalpha(:,:,1), dalpha(:,:,2))
    dC(:,:,2) = -dx*commutator(dalpha(:,:,1), 2*dalpha(:,:,3)+dC(:,:,1))/60

    dOmega = dalpha(:,:,1) + dalpha(:,:,3)/12 + &
            dx*commutator(-20*dalpha(:,:,1)-dalpha(:,:,3)+dC(:,:,1), dalpha(:,:,2)+dC(:,:,2))/240

    ! Finish

    return

  end subroutine eval_dOmega_GL6

!****

  function abscissa_magnus_GL2 (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(1)

    real(WP) :: dx

    ! Set up the abscissa for a 2nd-order Gauss-Legendre Magnus scheme

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    ! Finish

    return
    
  end function abscissa_magnus_GL2

!****

  function abscissa_magnus_GL4 (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(2)

    real(WP) :: dx

    ! Set up the abscissa for a 4th-order Gauss-Legendre Magnus scheme

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    ! Finish

    return

  end function abscissa_magnus_GL4

!****

  function abscissa_magnus_GL6 (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(3)

    real(WP) :: dx

    ! Set up the abscissa for a 6nd-order Gauss-Legendre Magnus scheme

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,0._WP,1._WP]*SQRT(15._WP)/10._WP)*dx

    ! Finish

    return

  end function abscissa_magnus_GL6

end module gyre_ivp_magnus
