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

  use gyre_bound
  use gyre_ext
  use gyre_match
  use gyre_status
  use gyre_sysmtx

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_t
     private
     class(${T}_bound_t), allocatable  :: bd_i
     class(${T}_bound_t), allocatable  :: bd_o
     type(${T}_bvp_seg_t), allocatable :: sg(:)
     type(${T}_match_t), allocatable   :: mt(:)
     class(${T}_sysmtx_t), allocatable :: sm
     real(WP)                          :: omega_min
     real(WP)                          :: omega_max
     integer                           :: n_e
     integer                           :: n_s
   contains 
     private
     procedure, public :: eval_discrim => eval_discrim_
     procedure         :: build_
     procedure         :: solve_
  end type ${T}_bvp_t

  ! Interfaces

  interface ${T}_bvp_t
     module procedure ${T}_bvp_t_
  end interface ${T}_bvp_t

  ! Access specifiers

  private

  public :: ${T}_bvp_t

contains

  function ${T}_bvp_t_ (bd_i, bd_o, sg, mt, sm, omega_min, omega_max) result (bp)

    class(${T}_bound_t), intent(in)   :: bd_i
    class(${T}_bound_t), intent(in)   :: bd_o
    class(${T}_bvp_seg_t), intent(in) :: sg
    class(${T}_match_t), intent(in)   :: mt
    class(${T}_sysmtx_t), intent(in)  :: sm
    real(WP), intent(in)              :: omega_min
    real(WP), intent(in)              :: omega_max
    type(${T}_bvp_t)                  :: bp

    integer :: n_e
    integer :: n_s

    ! Perform basic validations

    n_e = bd_i%n_e
    n_s = SIZE(sg)

    $CHECK_BOUNDS(bd_o%n_e,n_e)
    $CHECK_BOUNDS(bd_i%n+bd_o%n,n_e)

    $ASSERT_DEBUG(ALL(sg%n_e == n_e),Dimension mismatch)
    $ASSERT_DEBUG(ALL(mt%n_e == n_e),Dimension mismatch)

    $CHECK_BOUNDS(SIZE(mt),n_s-1)

    $ASSERT_DEBUG(SUM(sg%n-1)+n_s-1 == sm%n,Dimension mismatch)

    ! Construct the bvp_t

    allocate(bp%bd_i, SOURCE=bd_i)
    allocate(bp%bd_o, SOURCE=bd_o)

    allocate(bp%sg, SOURCE=sg)
    allocate(bp%mt, SOURCE=mt)

    allocate(bp%sm, SOURCE=sm)

    bp%n_e = n_e
    bp%n_s = SIZE(bp%sg)

    bp%omega_min = omega_min
    bp%omega_max = omega_max
    
    ! Finish
    
    return

  end function ${T}_bvp_t_

  !****

  subroutine eval_discrim_ (this, omega, discrim, status)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t), intent(out)    :: discrim
    integer, intent(out)             :: status

    ! Evaluate the discriminant as the determinant of the sysmtx

    if (REAL(omega, WP) >= this%omega_min .AND. REAL(omega, WP) <= this%omega_max) then

       call this%build_(omega)

       call this%sm%factorize()
       discrim = this%sm%det()

       status = STATUS_OK

    else

       status = STATUS_OMEGA_DOMAIN

    endif

    ! Finish

    return

  end subroutine eval_discrim_

  !****

  subroutine build_ (this, omega)

    class(${T}_bvp_t), target, intent(inout) :: this
    $TYPE(WP), intent(in)                    :: omega

    $TYPE(WP)                     :: E_l(this%n_e,this%n_e)
    $TYPE(WP)                     :: E_r(this%n_e,this%n_e)
    type(${T}_ext_t)              :: S_
    integer                       :: k
    integer                       :: s
    type(${T}_bvp_seg_t), pointer :: sg
    
    ! Build the sysmtx

    ! Boundary conditions

    call this%bd_i%build(omega, E_r, S_)
    call this%sm%set_B_i(E_r, S_)

    call this%bd_o%build(omega, E_l, S_)
    call this%sm%set_B_o(E_l, S_)

    ! Segments

    k = 1

    seg_loop : do s = 1, this%n_s

       sg => this%sg(s)

       !$OMP PARALLEL DO PRIVATE (E_l, E_r, S_) SCHEDULE (DYNAMIC)
       sub_loop : do i = 1, sg%n-1
          call sg%build(omega, i, E_l, E_r, S_)
          call this%sm%set_E(k+i-1, E_l, E_r, S_)
       end do sub_loop

       k = k + sg%n + 1

    end do seg_loop

    ! Match conditions

    k = 1

    match_loop : do s = 1, this%n_sg-1

       sg => this%sg(s)
       
       k = k + sg%n

       call this%mt(s)%build(omega, E_l, E_r, S_)
       call this%sm%set_E(k, E_l, E_r, S_)

       k = k + 1

    end do match_loop

    ! Finish

    return

  end subroutine build_

  !****

  subroutine solve_ (this, omega, discrim)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t), intent(out)    :: discrim

    $TYPE(WP) :: v(this%n_e*this%n)
    $TYPE(WP) :: y(this%n_e,this%n)

    ! Solve for the solution on the shooting grid, by calculating
    ! the null vector of the sysmtx

    call this%build_(omega)

    call this%sm%factorize()
    discrim = this%sm%det()
    v = this%sm%null_vec()

    y = RESHAPE(v, SHAPE(y))

    ! Copy y into the segments (this is probably inefficient; a
    ! better way would be to copy straight from v)

    k = 1

    seg_loop : do s = 1, this%n_s

       sg => this%sg(s)

       !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
       sub_loop : do i = 1, sg%n-1
          sg%y(:,i) = y(:,k+i-1)
       end do sub_loop

       k = k + sg%n + 1

    end do seg_loop

    ! Finish

    return

  end subroutine solve_

  ! !****

  ! subroutine recon_ (this, omega, x, x_ref, y, y_ref, discrim)

  !   class(${T}_bvp_t), intent(inout) :: this
  !   $TYPE(WP), intent(in)            :: omega
  !   real(WP), intent(in)             :: x(:)
  !   real(WP), intent(in)             :: x_ref
  !   $TYPE(WP), intent(out)           :: y(:,:)
  !   $TYPE(WP), intent(out)           :: y_ref(:)
  !   type(${T}_ext_t), intent(out)    :: discrim

  !   logical   :: same_grid
  !   integer   :: n
  !   integer   :: k
  !   integer   :: i_k(this%n)
  !   integer   :: i_a
  !   integer   :: i_b
  !   integer   :: k_ref
  !   $TYPE(WP) :: y_ref_(this%n_e,1)

  !   $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
  !   $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

  !   $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

  !   ! Reconstruct the solution
  !   ! Solve for y on the shooting grid

  !   call this%solve_(omega, discrim)

  !   ! Reconstruct the full solution in each subinterval

  !   if (SIZE(x) == SIZE(this%x)) then
  !      same_grid = ALL(x == this%x)
  !   else
  !      same_grid = .FALSE.
  !   endif

  !   if (same_grid) then

  !      y = this%y

  !   else

  !      n = SIZE(x)

  !      ! Determine where the shooting grid points lie in x

  !      call locate(x, this%x(1), i_k(1))

  !      do k = 2, this%n
  !         i_k(k) = i_k(k-1)
  !         call locate(x, this%x(k), i_k(k))
  !      end do

  !      ! Loop through subintervals

  !      !$OMP PARALLEL DO PRIVATE (i_a,i_b) SCHEDULE (DYNAMIC)
  !      sub_loop : do k = 1, this%n-1

  !         ! Select the points falling in the subinterval

  !         if (k == 1) then
  !            i_a = 1
  !            i_b = i_k(2)
  !         elseif (k == this%n-1) then
  !            i_a = i_k(this%n-1) + 1
  !            i_b = n
  !         else
  !            i_a = i_k(k) + 1
  !            i_b = i_k(k+1)
  !         endif

  !         ! Reconstruct in the subinterval

  !         if (i_b >= i_a) then
  !            call this%iv%recon(omega, this%x(k), this%x(k+1), this%y(:,k), this%y(:,k+1), &
  !                               x(i_a:i_b), y(:,i_a:i_b))
  !         endif

  !      end do sub_loop
          
  !   endif

  !   ! Reconstruct the solution at x_ref

  !   call locate(this%x, x_ref, k_ref)

  !   k_ref = MIN(MAX(k_ref, 1), this%n-1)
    
  !   call this%iv%recon(omega, this%x(k_ref), this%x(k_ref+1), this%y(:,k_ref), this%y(:,k_ref+1), &
  !                      [x_ref], y_ref_)

  !   y_ref = y_ref_(:,1)

  !   ! Finish

  !   return

  ! end subroutine recon_

end module gyre_${T}_bvp

$endsub
