! Incfile  : gyre_bvp
! Purpose  : boundary-value solver (template)
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_bvp

  ! Uses

  use core_kinds
  use core_order

  use gyre_bound
  use gyre_eqns
  use gyre_ext
  use gyre_ivp
  use gyre_model
  use gyre_rot
  use gyre_sysmtx

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_t
     class(model_t), pointer           :: ml => null()
     class(${T}_eqns_t), allocatable   :: eq
     class(${T}_bound_t), allocatable  :: bd
     class(${T}_ivp_t), allocatable    :: iv
     class(${T}_sysmtx_t), allocatable :: sm
     real(WP), allocatable             :: x(:)
     $TYPE(WP), allocatable            :: y(:,:)
     integer                           :: n
     integer                           :: n_e
   contains 
     private
     procedure, public :: discrim => discrim_
     procedure         :: build_
     procedure         :: solve_
     procedure, public :: recon => recon_
     procedure, public :: abscissa => abscissa_
  end type ${T}_bvp_t

  ! Interfaces

  interface ${T}_bvp_t
     module procedure ${T}_bvp_t_
  end interface ${T}_bvp_t

  ! Access specifiers

  private

  public :: ${T}_bvp_t

contains

  function ${T}_bvp_t_ (x, ml, eq, bd, iv, sm) result (bp)

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    class(${T}_eqns_t), intent(in)      :: eq
    class(${T}_bound_t), intent(in)     :: bd
    class(${T}_ivp_t), intent(in)       :: iv
    class(${T}_sysmtx_t), intent(in)    :: sm
    type(${T}_bvp_t), target            :: bp

    $CHECK_BOUNDS(bd%n_e,eq%n_e)

    $CHECK_BOUNDS(iv%n_e,eq%n_e)
    
    $CHECK_BOUNDS(sm%n_e,eq%n_e)
    $CHECK_BOUNDS(sm%n,SIZE(x)-1)

    ! Construct the bvp_t

    bp%n = SIZE(x)
    bp%n_e = eq%n_e

    bp%ml => ml
    
    allocate(bp%eq, SOURCE=eq)
    allocate(bp%bd, SOURCE=bd)
    allocate(bp%iv, SOURCE=iv)
    allocate(bp%sm, SOURCE=sm)

    bp%x = x
    allocate(bp%y(bp%n_e,bp%n))
    
    ! Finish
    
    return

  end function ${T}_bvp_t_

!****

  function discrim_ (this, omega) result (discrim)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t)                 :: discrim

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build_(omega)

    call this%sm%factorize()
    discrim = this%sm%det()

    ! Finish

    return

  end function discrim_

!****

  subroutine build_ (this, omega)

    class(${T}_bvp_t), target, intent(inout) :: this
    $TYPE(WP), intent(in)                    :: omega

    integer          :: k
    $TYPE(WP)        :: E_l(this%n_e,this%n_e)
    $TYPE(WP)        :: E_r(this%n_e,this%n_e)
    type(${T}_ext_t) :: S
    
    ! Build the sysmtx

    ! Boundary conditions

    call this%sm%set_B_i(this%bd%B_i(omega), ${T}_ext_t(1._WP))
    call this%sm%set_B_o(this%bd%B_o(omega), ${T}_ext_t(1._WP))

    ! Subintervals

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, S) SCHEDULE (DYNAMIC)
    sub_loop : do k = 1, this%n-1
       call this%iv%shoot(omega, this%x(k), this%x(k+1), E_l, E_r, S)
       call this%sm%set_E(k, E_l, E_r, S)
    end do sub_loop

    ! Finish

    return

  end subroutine build_

!****

  subroutine solve_ (this, omega, discrim)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t), intent(out)    :: discrim

    $TYPE(WP) :: v(this%n_e*this%n)

    ! Solve for the solution on the shooting grid, by calculating
    ! the null vector of the sysmtx

    call this%build_(omega)

    call this%sm%factorize()
    discrim = this%sm%det()
    v = this%sm%null_vec()

    this%y = RESHAPE(v, SHAPE(this%y))

    ! Finish

    return

  end subroutine solve_

!****

  subroutine recon_ (this, omega, x, x_ref, y, y_ref, discrim)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    real(WP), intent(in)             :: x(:)
    real(WP), intent(in)             :: x_ref
    $TYPE(WP), intent(out)           :: y(:,:)
    $TYPE(WP), intent(out)           :: y_ref(:)
    type(${T}_ext_t), intent(out)    :: discrim

    logical   :: same_grid
    integer   :: n
    integer   :: k
    integer   :: i_k(this%n)
    integer   :: i_a
    integer   :: i_b
    integer   :: k_ref
    $TYPE(WP) :: y_ref_(this%n_e,1)

    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

    ! Reconstruct the solution
    ! Solve for y on the shooting grid

    call this%solve_(omega, discrim)

    ! Reconstruct the full solution in each subinterval

    if (SIZE(x) == SIZE(this%x)) then
       same_grid = ALL(x == this%x)
    else
       same_grid = .FALSE.
    endif

    if (same_grid) then

       y = this%y

    else

       n = SIZE(x)

       ! Determine where the shooting grid points lie in x

       call locate(x, this%x(1), i_k(1))

       do k = 2, this%n
          i_k(k) = i_k(k-1)
          call locate(x, this%x(k), i_k(k))
       end do

       ! Loop through subintervals

       !$OMP PARALLEL DO PRIVATE (i_a,i_b) SCHEDULE (DYNAMIC)
       sub_loop : do k = 1, this%n-1

          ! Select the points falling in the subinterval

          if (k == 1) then
             i_a = 1
             i_b = i_k(2)
          elseif (k == this%n-1) then
             i_a = i_k(this%n-1) + 1
             i_b = n
          else
             i_a = i_k(k) + 1
             i_b = i_k(k+1)
          endif

          ! Reconstruct in the subinterval

          if (i_b >= i_a) then
             call this%iv%recon(omega, this%x(k), this%x(k+1), this%y(:,k), this%y(:,k+1), &
                                x(i_a:i_b), y(:,i_a:i_b))
          endif

       end do sub_loop
          
    endif

    ! Reconstruct the solution at x_ref

    call locate(this%x, x_ref, k_ref)

    k_ref = MIN(MAX(k_ref, 1), this%n-1)
    
    call this%iv%recon(omega, this%x(k_ref), this%x(k_ref+1), this%y(:,k_ref), this%y(:,k_ref+1), &
                       [x_ref], y_ref_)

    y_ref = y_ref_(:,1)

    ! Finish

    return

  end subroutine recon_

!****

  function abscissa_ (this) result (x)

    class(${T}_bvp_t), intent(in) :: this
    real(WP), allocatable         :: x(:)

    integer :: k
    integer :: n_cell(this%n-1)
    integer :: i

    ! Determine the abscissa used for shooting

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
    count_loop : do k = 1, this%n-1
       n_cell(k) = SIZE(this%iv%abscissa(this%x(k), this%x(k+1)))
    end do count_loop

    allocate(x(SUM(n_cell)))

    i = 1

    cell_loop : do k = 1, this%n-1
       x(i:i+n_cell(k)-1) = this%iv%abscissa(this%x(k), this%x(k+1))
       i = i + n_cell(k)
    end do cell_loop

    $CHECK_BOUNDS(i,SIZE(x)+1)

    ! Finish

    return

  end function abscissa_

end module gyre_${T}_bvp

$endsub
