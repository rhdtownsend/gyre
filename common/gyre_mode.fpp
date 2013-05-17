! Module   : gyre_mode
! Purpose  : mode data
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

module gyre_mode

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

  type :: mode_t
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
     procedure, public :: delrho
     procedure, public :: delT
     procedure, public :: dE_dx
     procedure, public :: dW_dx
     procedure, public :: C
     procedure, public :: E
     procedure, public :: E_norm
     procedure, public :: W
     procedure, public :: omega_im
  end type mode_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_ef
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: mode_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, bc, tc, op, omega, x, y)

    class(mode_t), intent(out)                     :: this
    class(base_coeffs_t), intent(in)               :: bc
    class(therm_coeffs_t), intent(in), allocatable :: tc
    type(oscpar_t), intent(in)                     :: op
    complex(WP), intent(in)                        :: omega
    real(WP), intent(in)                           :: x(:)
    complex(WP), intent(in)                        :: y(:,:)

    real(WP) :: phase

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Initialize the mode

    allocate(this%bc, SOURCE=bc)
    if(ALLOCATED(tc)) allocate(this%tc, SOURCE=tc)

    this%op = op

    this%x = x
    this%y = y

    this%omega = omega

    this%n = SIZE(this%x)

    ! Normalize by the mode inertia, and so that y(1,n) is real

    phase = ATAN2(AIMAG(this%y(1,this%n)), REAL(this%y(1,this%n)))

    this%y = this%y/SQRT(this%E())*EXP(CMPLX(0._WP, -phase, KIND=WP))

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_ef (this, root_rank)

    class(mode_t), intent(inout) :: this
    integer, intent(in)          :: root_rank

    ! Broadcast the mode

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

    class(mode_t), intent(in) :: this
    integer, intent(out)      :: n_p
    integer, intent(out)      :: n_g

    real(WP) :: y_1(this%n)
    real(WP) :: y_2(this%n)
    logical  :: inner_ext
    integer  :: i
    real(WP) :: y_2_cross

    ! Classify the mode using the Eckart-Scuflaire-Osaki scheme

    y_1 = REAL(this%y(1,:))
    y_2 = REAL(this%y(2,:))

    if(this%op%l == 0) then
       n_p = 1
       n_g = 0
       inner_ext = .FALSE.
    else
       n_p = 0
       n_g = 0
    endif
 
    x_loop : do i = 1,this%n-1

       ! If this is a radial mode, and extremum in y_1 hasn't yet been
       ! reached, skip (this is to deal with noisy near-zero solutions
       ! at the origin)

       if(this%op%l == 0) then
          if(i > 1 .AND. .NOT. inner_ext) inner_ext = ABS(y_1(i)) > ABS(y_1(i-1)) .AND. ABS(y_1(i)) > ABS(y_1(i+1))
          if(.NOT. inner_ext) cycle x_loop
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

    class(mode_t), intent(in) :: this
    complex(WP)               :: xi_r(this%n)
    
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

    class(mode_t), intent(in) :: this
    complex(WP)               :: xi_h(this%n)
    
    ! Calculate the radial displacement perturbation in units of
    ! R_star

    associate (c_1 => this%bc%c_1(this%x), l => this%op%l, omega => this%omega)

      if(l /= 0) then

         if(l /= 1) then

            where (this%x > 0._WP)
               xi_h = this%y(2,:)*this%x**(l-1)/(c_1*omega**2)
            elsewhere
               xi_h = 0._WP
            end where

         else

            xi_h = this%y(2,:)/(c_1*omega**2)

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

    class(mode_t), intent(in) :: this
    complex(WP)               :: phip(this%n)
    
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

    class(mode_t), intent(in) :: this
    complex(WP)               :: dphip_dx(this%n)
    
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

    class(mode_t), intent(in) :: this
    complex(WP)               :: delS(this%n)

    ! Calculate the Lagrangian specific entropy perturbation in units
    ! of c_p

    associate(l => this%op%l)

      where (this%x /= 0._WP)
         delS = this%y(5,:)*this%x**(l-2)
      elsewhere
         delS = 0._WP
      end where

    end associate
         
    ! Finish

    return

  end function delS

!****

  function delL (this)

    class(mode_t), intent(in)  :: this
    complex(WP)                :: delL(this%n)

    ! Calculate the Lagrangian luminosity perturbation in units of R_star

    associate (l => this%op%l)

      delL = this%y(6,:)*this%x**(l+1)

    end associate

    ! Finish

    return

  end function delL

!****

  function delp (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delp(this%n)

    ! Calculate the Lagrangian pressure perturbation in units of p

    associate (V => this%bc%V(this%x), pi_c => this%bc%pi_c(), l => this%op%l)

      if(l > 0) then

         where (this%x /= 0._WP)
            delp = V*(this%y(2,:) - this%y(1,:) - this%y(3,:))*this%x**(l-2)
         elsewhere
            delp = 0._WP
         endwhere

      else

         where (this%x /= 0._WP)
            delp = V*(this%y(2,:) - this%y(1,:) - this%y(3,:))*this%x**(l-2)
         elsewhere
            delp = pi_c*(this%y(2,:) - this%y(1,:) - this%y(3,:))
         endwhere

      endif

    end associate

    ! Finish

    return

  end function delp

!****

  function delrho (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delrho(this%n)

    ! Calculate the Lagrangian density perturbation in units of rho

    associate (Gamma_1 => this%bc%Gamma_1(this%x))

      delrho = this%delp()/Gamma_1 - this%delS()

    end associate

    ! Finish

    return

  end function delrho

!****

  function delT (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delT(this%n)

    ! Calculate the Lagrangian temperature perturbation in units of T

    associate (nabla_ad => this%bc%nabla_ad(this%x))

      delT = nabla_ad*this%delp() + this%delS()

    end associate

    ! Finish

    return

  end function delT

!****

  function dE_dx (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: dE_dx(this%n)

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)
    
    ! Calculate the differential mode inertia (Aerts et al. 2010,
    ! eqn. 3.139) in units of M_star R_star**2

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(U => this%bc%U(this%x), c_1 => this%bc%c_1(this%x), &
              l => this%op%l)
      dE_dx = 4._WP*PI*(ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*this%x**2/c_1
    end associate

    ! Finish

    return

  end function dE_dx

!****

  function dW_dx (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: dW_dx(this%n)

    $ASSERT(ALLOCATED(this%tc),No therm_coeffs data)

    ! Calculate the differential work in units of G M_star^2/R_star
    ! t_dyn/t_KH = t_dyn L_*.  The entropy-based expression for the
    ! work is used (cf. Unno et al.  1989, eqn. 25.9); the additional
    ! factor of 4 pi in the denominator comes from averaging over
    ! solid angle

    associate(c_thm => this%tc%c_thm(this%x))

      dW_dx = -PI*AIMAG(CONJG(this%delT())*this%delS())*c_thm*this%x**2/(4._WP*PI)

    end associate

    ! Finish

    return

  end function dW_dx

!****

  function C (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: C

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)

    ! Calculate the first-order rotational splitting coefficient (Unno
    ! et al. 1989, eqn. 19.46)

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(x => this%x, U => this%bc%U(this%x), c_1 => this%bc%c_1(this%x), &
              l => this%op%l)

      C = integrate(this%x, (2._WP*REAL(xi_r*CONJG(xi_h)) + ABS(xi_h)**2)*U*x**2/c_1)/ &
          integrate(this%x, (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*x**2/c_1)

    end associate

    ! Finish

    return

  end function C

!****

  function E (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E
    
    ! Calculate the total mode inertia (Aerts et al. 2010, eqn. 3.139)
    ! in units of M_star R_star**2

    E = integrate(this%x, this%dE_dx())

    ! Finish

    return

  end function E

!****

  function E_norm (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)
    real(WP)    :: A2

    ! Calculate the normalized mode inertia (Aerts et al. 2010,
    ! eqn. 3.140)

    E = this%E()

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(l => this%op%l)

      A2 = ABS(xi_r(this%n))**2 + l*(l+1)*ABS(xi_h(this%n))**2

      if(A2 == 0._WP) then
         $WARN(Surface amplitude is zero; not normalizing inertia)
         E_norm = E
      else
         E_norm = E/A2
      endif

    end associate

    ! Finish

    return

  end function E_norm

!*****

  function W (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: W
    
    ! Calculate the total work

    W = integrate(this%x, this%dW_dx())

    ! Finish

    return

  end function W

!****

  function omega_im (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: omega_im

    integer  :: i_trans
    integer  :: i
    real(WP) :: dW_dx(this%n)
    real(WP) :: W
    
    ! Estimate the imaginary part of omega by integrating the work
    ! function out to the thermal transition point

    ! First locate the point

    i_trans = this%n

    do i = this%n-1,1,-1

       associate(tau_thm => this%tc%tau_thm(this%x(i)))

         if(REAL(this%omega)*tau_thm/TWOPI > 1._WP) then
            i_trans = i
            exit
         endif

       end associate

    enddo

    i_trans = this%n

    ! Do the integration

    dW_dx = this%dW_dx()

    W = integrate(this%x(:i_trans), dW_dx(:i_trans))

    ! Calculate omega_im

    omega_im = -REAL(this%omega)*W/this%E()

    ! Finish

    return

  end function omega_im

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

end module gyre_mode
