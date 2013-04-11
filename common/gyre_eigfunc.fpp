! Module   : gyre_eigfunc
! Purpose  : eigenfunction data
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

module gyre_eigfunc

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup

  use gyre_bvp
  use gyre_base_coeffs
  use gyre_therm_coeffs
  $if($MPI)
  use gyre_base_coeffs_mpi
  use gyre_therm_coeffs_mpi
  $endif
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: eigfunc_t
     class(base_coeffs_t), allocatable  :: bc
     class(therm_coeffs_t), allocatable :: tc
     type(oscpar_t)                     :: op
     real(WP), allocatable              :: x(:)
     complex(WP), allocatable           :: y(:,:)
     complex(WP)                        :: omega
     integer                            :: n
   contains
     private
     procedure, public :: init
     procedure, public :: classify
     procedure, public :: xi_r
     procedure, public :: xi_h
     procedure, public :: phip
     procedure, public :: dphip_dx
     procedure, public :: delS
     procedure, public :: delL
     procedure, public :: delp
     procedure, public :: delT
     procedure, public :: dK_dx
     procedure, public :: dW_dx
     procedure, public :: E
     procedure, public :: K
     procedure, public :: W
  end type eigfunc_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_ef
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: eigfunc_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, bc, tc, op, omega, x, y)

    class(eigfunc_t), intent(out)                  :: this
    class(base_coeffs_t), intent(in)               :: bc
    class(therm_coeffs_t), intent(in), allocatable :: tc
    type(oscpar_t), intent(in)                     :: op
    complex(WP), intent(in)                        :: omega
    real(WP), intent(in)                           :: x(:)
    complex(WP), intent(in)                        :: y(:,:)

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Initialize the eigfunc

    allocate(this%bc, SOURCE=bc)
    if(ALLOCATED(tc)) allocate(this%tc, SOURCE=tc)

    this%op = op

    this%x = x
    this%y = y

    this%omega = omega

    this%n = SIZE(this%x)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_ef (this, root_rank)

    class(eigfunc_t), intent(inout) :: this
    integer, intent(in)             :: root_rank

    ! Broadcast the eigfunc

    call bcast_alloc(this%bc, root_rank)
    call bcast_alloc(this%tc, root_rank)

    call bcast(this%op, root_rank)

    call bcast_alloc(this%x, root_rank)
    call bcast_alloc(this%y, root_rank)

    call bcast(this%omega, root_rank)

    call bcast(this%n, root_rank)

  end subroutine bcast_ef

  $endif

!****

  subroutine classify (this, n_p, n_g)

    class(eigfunc_t), intent(in) :: this
    integer, intent(out)         :: n_p
    integer, intent(out)         :: n_g

    real(WP) :: y_1(this%n)
    real(WP) :: y_2(this%n)
    logical  :: inner_ext
    integer  :: i
    real(WP) :: y_2_cross

    ! Classify the eigenfunction using the Cowling-Scuflaire scheme

    y_1 = REAL(this%y(1,:))
    y_2 = REAL(this%y(2,:))

    n_p = 0
    n_g = 0
 
    inner_ext = ABS(y_1(1)) > ABS(y_1(2))

    x_loop : do i = 2,this%n-1

       ! If the innermost extremum in y_1 hasn't yet been reached,
       ! skip

       if(.NOT. inner_ext) then
          inner_ext = ABS(y_1(i)) > ABS(y_1(i-1)) .AND. ABS(y_1(i)) > ABS(y_1(i+1))
          cycle x_loop
       endif

       ! Look for a node in xi_r

       if(y_1(i) >= 0._WP .AND. y_1(i+1) < 0._WP) then

          y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

          if(y_2_cross >= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       elseif(y_1(i) <= 0._WP .AND. y_1(i+1) > 0._WP) then

         y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

          if(y_2_cross <= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       endif

    end do x_loop

    ! Finish

    return

  end subroutine classify

!****

  function xi_r (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: xi_r(this%n)
    
    ! Calculate the radial displacement perturbation in units of
    ! R_star

    associate (l => this%op%l)

      if(l /= 1) then

         where (this%x > 0._WP)
            xi_r = this%y(1,:)*this%x**(l-1)
         elsewhere
            xi_r = 0._WP
         end where

      else

         xi_r = this%y(1,:)

      end if

    end associate

    ! Finish

    return

  end function xi_r

!****

  function xi_h (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: xi_h(this%n)
    
    ! Calculate the radial displacement perturbation in units of
    ! R_star

    associate (c_1 => this%bc%c_1(this%x), l => this%op%l)

      if(l /= 0) then

         if(l /= 1) then

            where (this%x > 0._WP)
               xi_h = this%y(2,:)*this%x**(l-1)/(c_1*this%omega**2)
            elsewhere
               xi_h = 0._WP
            end where

         else

            xi_h = this%y(2,:)/(c_1*this%omega**2)

         endif

      else

         xi_h = 0._WP

      end if

    end associate

    ! Finish

    return

  end function xi_h

!****

  function phip (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: phip(this%n)
    
    ! Calculate the Eulerian gravitational potential perturbation in units of
    ! G M_star / R_star

    associate (c_1 => this%bc%c_1(this%x), l => this%op%l)

      phip = this%y(3,:)*this%x**l/c_1

    end associate

    ! Finish

    return

  end function phip

!****

  function dphip_dx (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: dphip_dx(this%n)
    
    ! Calculate the Eulerian gravity perturbation in units of G M_star
    ! / R_star

    associate (c_1 => this%bc%c_1(this%x), l => this%op%l)

      if(l /= 1) then

         where (this%x > 0._WP)
            dphip_dx = this%y(4,:)*this%x**(l-1)/c_1
         elsewhere
            dphip_dx = 0._WP
         end where

      else

         dphip_dx = this%y(4,:)/c_1

      end if

    end associate

    ! Finish

    return

  end function dphip_dx

!****

  function delS (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: delS(this%n)
    
    ! Calculate the Lagrangian specific entropy perturbation in units
    ! of c_p

    associate (l => this%op%l)

      delS = this%y(5,:)*this%x**l

    end associate

    ! Finish

    return

  end function delS

!****

  function delL (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: delL(this%n)
    
    ! Calculate the Lagrangian luminosity perturbation in units of R_star

    associate (l => this%op%l)

      delL = this%y(6,:)*this%x**(l+1)

    end associate

    ! Finish

    return

  end function delL

!****

  function delp (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: delp(this%n)

    ! Calculate the Lagrangian pressure perturbation in units of p

    associate (V_x2 => this%bc%V_x2(this%x), l => this%op%l)

      delp = V_x2*(this%y(2,:) - this%y(1,:) - this%y(3,:))*this%x**l

    end associate

    ! Finish

    return

  end function delp

!****

  function delT (this)

    class(eigfunc_t), intent(in) :: this
    complex(WP)                  :: delT(this%n)

    ! Calculate the Lagrangian temperature perturbation in units of T

    associate (nabla_ad => this%bc%nabla_ad(this%x))

      delT = nabla_ad*this%delp() + this%delS()

    end associate

    ! Finish

    return

  end function delT

!****

  function dK_dx (this)

    class(eigfunc_t), intent(in) :: this
    real(WP)                     :: dK_dx(this%n)

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)
    
    ! Calculate the differential kinetic energy in units of G M_star^2/R_star

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(U => this%bc%U(this%x), c_1 => this%bc%c_1(this%x), &
              l => this%op%l)
      dK_dx = (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h))*U*this%x**2/c_1
    end associate

    ! Finish

    return

  end function dK_dx

!****

  function dW_dx (this)

    class(eigfunc_t), intent(in) :: this
    real(WP)                     :: dW_dx(this%n)

    $ASSERT(ALLOCATED(this%tc),No therm_coeffs data)

    ! Calculate the differential work in units of G M_star^2/R_star
    ! t_dyn/t_KH.  The entropy-based expression for the work is used
    ! (cf. Unno et al.  1989, eqn. 25.9); the additional factor of 4
    ! pi in the denominator comes from averaging over solid angle

    associate(c_thm => this%tc%c_thm(this%x))

      dW_dx = -PI*AIMAG(CONJG(this%delT())*this%delS())*c_thm*this%x**2/(4._WP*PI)

    end associate

    ! Finish

    return

  end function dW_dx

!*****

  function E (this)

    class(eigfunc_t), intent(in) :: this
    real(WP)                     :: E

    real(WP)    :: K
    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)
    real(WP)    :: A2

    ! Calculate the normalized mode inertia, using the expression
    ! given by Christensen-Dalsgaard (2011, arXiv:1106.5946, his
    ! eqn. 2.32)

    K = this%K()

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(l => this%op%l)

      A2 = ABS(xi_r(this%n))**2 + l*(l+1)*ABS(xi_h(this%n))**2

      if(A2 == 0._WP) then
         $WARN(Surface amplitude is zero; not normalizing inertia)
         E = K
      else
         E = K/A2
      endif

    end associate

    ! Finish

    return

  end function E

!*****

  function K (this)

    class(eigfunc_t), intent(in) :: this
    real(WP)                     :: K
    
    ! Calculate the kinetic energy

    K = integrate(this%x, this%dK_dx())

    ! Finish

    return

  end function K

!*****

  function W (this)

    class(eigfunc_t), intent(in) :: this
    real(WP)                     :: W
    
    ! Calculate the total work

    W = integrate(this%x, this%dW_dx())

    ! Finish

    return

  end function W

!****

  function integrate (x, y) result (int_y)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: int_y

    integer :: n

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Integrate y(x) using trapezoidal quadrature

    n = SIZE(x)

    int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))

    ! Finish

    return

  end function integrate

end module gyre_eigfunc
