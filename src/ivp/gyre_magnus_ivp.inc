! Incfile  : gyre_magnus_ivp
! Purpose  : initial-value solvers (Magnus method, template)
!
! Copyright 2013-2014 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_magnus_ivp

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_ext
  use gyre_ivp
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GL2_SCHEME = 1
  integer, parameter :: GL4_SCHEME = 2
  integer, parameter :: GL6_SCHEME = 3

  ! Derived-type definitions

  type, extends (${T}_ivp_t) :: ${T}_magnus_ivp_t
     private
     class(${T}_eqns_t), allocatable :: eq
     integer                         :: scheme
   contains
     private
     procedure, public :: shoot => shoot_
     procedure, public :: recon => recon_
     procedure, public :: dOmega => dOmega_
  end type ${T}_magnus_ivp_t

  ! Interfaces

  interface ${T}_magnus_ivp_t
     module procedure ${T}_magnus_ivp_t_
  end interface ${T}_magnus_ivp_t

  ! Access specifiers

  private

  public :: ${T}_magnus_ivp_t

  ! Procedures

contains

  function ${T}_magnus_ivp_t_ (eq, scheme) result (iv)

    class(${T}_eqns_t), intent(in) :: eq
    character(*), intent(in)       :: scheme
    type(${T}_magnus_ivp_t)        :: iv

    ! Construct the magnus_ivp_t

    allocate(iv%eq, SOURCE=eq)

    select case (scheme)
    case ('GL2')
       iv%scheme = GL2_SCHEME
    case ('GL4')
       iv%scheme = GL4_SCHEME
    case ('GL6')
       iv%scheme = GL6_SCHEME
    case default
       $ABORT(Invalid scheme)
    end select

    iv%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_magnus_ivp_t_

!****

  subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: S

    logical, parameter :: UPWIND = .TRUE.
    
    $TYPE(WP)   :: dOmega(this%n_e,this%n_e)
    complex(WP) :: lambda(this%n_e)
    complex(WP) :: V_l(this%n_e,this%n_e)
    complex(WP) :: V_r(this%n_e,this%n_e)
    real(WP)    :: dx
    integer     :: i
    complex(WP) :: V_pos(this%n_e,this%n_e)
    complex(WP) :: V_neg(this%n_e,this%n_e)
    complex(WP) :: E_l_(this%n_e,this%n_e)
    complex(WP) :: E_r_(this%n_e,this%n_e)

    ! Evaluate the Magnus slope matrix

    dOmega = this%dOmega(omega, x_a, x_b)

    ! Decompose it

    call eigen_decompose(dOmega, lambda, V_l, V_r)

    ! Set up the shooting matrices and scales

    dx = x_b - x_a

    $block

    $if ($DOUBLE_PRECISION)
    $local $X Z
    $else
    $local $X C
    $endif

    if (UPWIND) then

       ! Use 'upwinding' for stability

       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_pos(1,i), 1)
          if (REAL(lambda(i)) >= 0._WP) then
             call ${X}SCAL(this%n_e, EXP(-lambda(i)*dx), V_pos(1,i), 1)
          endif
       end do
    
       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          if (REAL(lambda(i)) < 0._WP) then
             call ${X}SCAL(this%n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
          endif
       end do

       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(-1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l_, this%n_e)

       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_pos, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_r_, this%n_e)

       S = exp(${T}_ext_t(SUM(lambda, MASK=REAL(lambda) >= 0._WP)*dx))

    else

       ! Use the [TowTei2013] expression (their eqns. 4, 16, 18)

       do i = 1,this%n_e
          call ${X}COPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          call ${X}SCAL(this%n_e, EXP(lambda(i)*dx), V_neg(1,i), 1)
       end do
    
       call ${X}GEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l_, this%n_e)

       do i = 1, this%n_e
          E_r_(:,i) = 0._WP
          E_r_(i,i) = -1._WP
       end do

       S = ${T}_ext_t(1._WP)

    endif

    $endblock

    E_l = E_l_
    E_r = E_r_

    ! Finish

    return

  end subroutine shoot_

!****

  subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(in)                :: y_a(:)
    $TYPE(WP), intent(in)                :: y_b(:)
    real(WP), intent(in)                 :: x(:)
    $TYPE(WP), intent(out)               :: y(:,:)

    $TYPE(WP)   :: dOmega(this%n_e,this%n_e)
    complex(WP) :: lambda(this%n_e)
    complex(WP) :: V_l(this%n_e,this%n_e)
    complex(WP) :: V_r(this%n_e,this%n_e)
    integer     :: i
    complex(WP) :: exp_a(this%n_e)
    complex(WP) :: exp_b(this%n_e)

    ! Evaluate the Magnus slope matrix

    dOmega = this%dOmega(omega, x_a, x_b)

    ! Decompose it

    call eigen_decompose(dOmega, lambda, V_l, V_r)

    ! Do the stabilized (both-boundaries) Magnus reconstruction

    recon_loop : do i = 1,SIZE(x)

       exp_a = MERGE(EXP(lambda*(x(i) - x_a)), CMPLX(0._WP, KIND=WP), REAL(lambda) < 0._WP .EQV. x_b > x_a)
       exp_b = MERGE(EXP(lambda*(x(i) - x_b)), CMPLX(0._WP, KIND=WP), REAL(lambda) >= 0._WP .EQV. x_b > x_a)
       
       y(:,i) = MATMUL(V_r, MATMUL(diagonal_matrix(exp_a), MATMUL(V_l, y_a)) + &
                            MATMUL(diagonal_matrix(exp_b), MATMUL(V_l, y_b)))

    end do recon_loop

    ! Finish

    return

  end subroutine recon_

!****

  function dOmega_ (this, omega, x_a, x_b) result (dOmega)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP)                            :: dOmega(this%n_e,this%n_e)

    ! Evaluate the Magnus slope matrix

    select case (this%scheme)
    case (GL2_SCHEME)
       dOmega = dOmega_GL2_(this, omega, x_a, x_b)
    case (GL4_SCHEME)
       dOmega = dOmega_GL4_(this, omega, x_a, x_b)
    case (GL6_SCHEME)
       dOmega = dOmega_GL6_(this, omega, x_a, x_b)
    case default
       $ABORT(Invalid scheme)
    end select

    ! Finish

    return

  end function dOmega_

!****

  function dOmega_GL2_ (this, omega, x_a, x_b) result (dOmega)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP)                            :: dOmega(this%n_e,this%n_e)

    real(WP)  :: dx
    real(WP)  :: x(1)
    $TYPE(WP) :: A(this%n_e,this%n_e,1)

    ! Evaluate the GL2 Magnus slope matrix

    ! Calculate the RHS matrix

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    A(:,:,1) = this%eq%A(x(1), omega)

    ! Set up the slope matrix

    dOmega = A(:,:,1)

    ! Finish

    return

  end function dOmega_GL2_

!****

  function dOmega_GL4_ (this, omega, x_a, x_b) result (dOmega)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP)                            :: dOmega(this%n_e,this%n_e)

    real(WP)  :: dx
    real(WP)  :: x(2)
    $TYPE(WP) :: A(this%n_e,this%n_e,2)
    $TYPE(WP) :: dalpha(this%n_e,this%n_e,2)

    ! Evaluate the GL4 Magnus slope matrix

    ! Calculate the RHS matrices

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    A(:,:,1) = this%eq%A(x(1), omega)
    A(:,:,2) = this%eq%A(x(2), omega)

    ! Set up the Magnus slope matrix (Blanes et al. 2009, eqns. 243
    ! and 253; note that the 12 in the denominator of their expression
    ! for alpha_2 is erroneous)

    dalpha(:,:,1) = 0.5_WP*(A(:,:,1) + A(:,:,2))
    dalpha(:,:,2) = SQRT(3._WP)*(A(:,:,2) - A(:,:,1))

    dOmega = dalpha(:,:,1) - dx*commutator(dalpha(:,:,1), dalpha(:,:,2))/12._WP
    
    ! Finish

    return

  end function dOmega_GL4_

!****

  function dOmega_GL6_ (this, omega, x_a, x_b) result (dOmega)

    class(${T}_magnus_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP)                            :: dOmega(this%n_e,this%n_e)

    real(WP)  :: dx
    real(WP)  :: x(3)
    $TYPE(WP) :: A(this%n_e,this%n_e,3)
    $TYPE(WP) :: dalpha(this%n_e,this%n_e,3)
    $TYPE(WP) :: dC(this%n_e,this%n_e,2)

    ! Evaluate the GL6 Magnus slope matrix

    ! Calculate the RHS matrices

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,0._WP,1._WP]*SQRT(15._WP)/10._WP)*dx

    A(:,:,1) = this%eq%A(x(1), omega)
    A(:,:,2) = this%eq%A(x(2), omega)
    A(:,:,3) = this%eq%A(x(3), omega)

    ! Set up the Magnus slope matrix (Blanes et al. 2009, eqns. 251
    ! and 257)

    dalpha(:,:,1) = A(:,:,2)
    dalpha(:,:,2) = SQRT(15._WP)*(A(:,:,3) - A(:,:,1))/3
    dalpha(:,:,3) = 10*(A(:,:,3) - 2*A(:,:,2) + A(:,:,1))/3

    dC(:,:,1) = dx*commutator(dalpha(:,:,1), dalpha(:,:,2))
    dC(:,:,2) = -dx*commutator(dalpha(:,:,1), 2*dalpha(:,:,3)+dC(:,:,1))/60

    dOmega = dalpha(:,:,1) + dalpha(:,:,3)/12 + &
            dx*commutator(-20*dalpha(:,:,1)-dalpha(:,:,3)+dC(:,:,1), dalpha(:,:,2)+dC(:,:,2))/240

    ! Finish

    return

  end function dOmega_GL6_

end module gyre_${T}_magnus_ivp

$endsub
