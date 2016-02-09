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
  use gyre_num_par
  use gyre_status
  use gyre_sysmtx
  use gyre_sysmtx_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_t
     private
     class(${T}_bound_t), allocatable          :: bd_i
     class(${T}_bound_t), allocatable          :: bd_o
     type(${T}_bvp_seg_t), allocatable, public :: bs(:)
     type(${T}_match_t), allocatable           :: mt(:)
     class(${T}_sysmtx_t), allocatable         :: sm
     real(WP)                                  :: omega_min
     real(WP)                                  :: omega_max
     integer, public                           :: n_e
     integer, public                           :: n_s
   contains 
     private
     procedure, public :: eval_discrim => eval_discrim_
     procedure         :: build_
     procedure, public :: recon => recon_
  end type ${T}_bvp_t

  ! Interfaces

  interface ${T}_bvp_t
     module procedure ${T}_bvp_t_
  end interface ${T}_bvp_t

  ! Access specifiers

  private

  public :: ${T}_bvp_t

contains

  function ${T}_bvp_t_ (bd_i, bd_o, bs, mt, nm_p, omega_min, omega_max) result (bp)

    class(${T}_bound_t), intent(in)   :: bd_i
    class(${T}_bound_t), intent(in)   :: bd_o
    class(${T}_bvp_seg_t), intent(in) :: bs
    class(${T}_match_t), intent(in)   :: mt
    type(num_par_t), intent(in)       :: nm_p
    real(WP), intent(in)              :: omega_min
    real(WP), intent(in)              :: omega_max
    type(${T}_bvp_t)                  :: bp

    integer :: n_e
    integer :: n_s
    integer :: n_i
    integer :: n_o

    ! Perform basic validations

    n_e = bd_i%n_e
    n_s = SIZE(bs)

    n_i = bd_i%n
    n_o = bd_o%n

    $CHECK_BOUNDS(bd_o%n_e,n_e)
    $CHECK_BOUNDS(bd_i%n+bd_o%n,n_e)

    $ASSERT_DEBUG(ALL(bs%n_e == n_e),Dimension mismatch)
    $ASSERT_DEBUG(ALL(mt%n_e == n_e),Dimension mismatch)

    $CHECK_BOUNDS(SIZE(mt),n_s-1)

    ! Construct the bvp_t

    allocate(bp%bd_i, SOURCE=bd_i)
    allocate(bp%bd_o, SOURCE=bd_o)

    allocate(bp%bs, SOURCE=bs)
    allocate(bp%mt, SOURCE=mt)

    allocate(bp%sm, SOURCE=${T}_sysmtx_t(SUM(bs%n)-1, n_e, n_i, n_o, nm_p)

    bp%n_e = n_e
    bp%n_s = n_s

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

    $TYPE(WP)            :: E_l(this%n_e,this%n_e)
    $TYPE(WP)            :: E_r(this%n_e,this%n_e)
    type(${T}_ext_t)     :: S_
    integer              :: k
    integer              :: s
    integer, allocatable :: i_seg(:)
    integer              :: i_dbl(this%pt%n_dbl)
    
    ! Build the sysmtx

    ! Boundary conditions

    call this%bd_i%build(omega, E_r, S_)
    call this%sm%set_B_i(E_r, S_)

    call this%bd_o%build(omega, E_l, S_)
    call this%sm%set_B_o(E_l, S_)

    ! Segments

    ! WHERE IS THE PARTITIONER SET UP?

    seg_loop : do s = 1, this%n_s

       i_seg = this%pt%i_seg(s)

       associate (bs => this%bs(s))

         !$OMP PARALLEL DO PRIVATE (E_l, E_r, S_) SCHEDULE (DYNAMIC)
         sub_loop : do i = 1, bs%n-1
            call bs%build(omega, i, E_l, E_r, S_)
            call this%sm%set_E(i_seg(i), E_l, E_r, S_)
         end do sub_loop

       end associate

    end do seg_loop

    ! Match conditions

    i_dbl = this%pt%i_dbl()

    match_loop : do s = 1, this%n_s-1

       associate (bs => this%bs(s))

         call this%mt(s)%build(omega, E_l, E_r, S_)
         call this%sm%set_E(i_dbl(s), E_l, E_r, S_)

       end associate

    end do match_loop

    ! Finish

    return

  end subroutine build_

  !****

  subroutine recon_ (this, omega, discrim)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t), intent(out)    :: discrim

    $TYPE(WP) :: v(this%n_e*this%n)
    $TYPE(WP) :: y(this%n_e,this%n)

    ! Reconstruct the solution, by calculating the null vector of the
    ! sysmtx

    call this%build_(omega)

    call this%sm%factorize()
    discrim = this%sm%det()
    v = this%sm%null_vec()

    y = RESHAPE(v, SHAPE(y))

    ! Copy y into the segments

    DONT YET KNOW HOW TO COPY INTO SEGMENTS

    ! Finish

    return

  end subroutine recon_

end module gyre_${T}_bvp

$endsub
