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

  use gyre_coeffs
  $if($MPI)
  use gyre_coeffs_mpi
  $endif
  use gyre_ext_arith
  use gyre_oscpar
  use gyre_grid
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_t
     class(coeffs_t), pointer :: cf => null()
     type(oscpar_t)           :: op
     type(ext_real_t)         :: chi
     real(WP), allocatable    :: x(:)
     complex(WP), allocatable :: y(:,:)
     real(WP)                 :: x_ref
     complex(WP)              :: y_ref(6)
     complex(WP)              :: omega
     integer                  :: n_pg
     integer                  :: n_p
     integer                  :: n_g
     integer                  :: n
     integer                  :: n_iter
   contains
     private
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $endif
     procedure, public :: freq
     procedure, public :: xi_r
     procedure, public :: xi_h
     procedure, public :: Yt_1
     procedure, public :: Yt_2
     procedure, public :: phip
     procedure, public :: dphip_dx
     procedure, public :: delS
     procedure, public :: delS_en
     procedure, public :: delL
     procedure, public :: delL_rd
     procedure, public :: delp
     procedure, public :: delrho
     procedure, public :: delT
     procedure, public :: dE_dx
     procedure, public :: dW_dx
     procedure, public :: I_0
     procedure, public :: I_1
     procedure, public :: prop_type
     procedure, public :: K
     procedure, public :: beta
     procedure, public :: E
     procedure, public :: E_norm
     procedure, public :: W
     procedure, public :: omega_im
     procedure, public :: xi_r_ref
     procedure, public :: xi_h_ref
  end type mode_t

  ! Interfaces

  interface mode_t
     module procedure init_md
  end interface mode_t

  $if ($MPI)
  interface bcast
     module procedure bcast_md
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

  function init_md (cf, op, omega, x, y, x_ref, y_ref, chi, n_iter) result (md)

    class(coeffs_t), pointer, intent(in) :: cf
    type(oscpar_t), intent(in)           :: op
    complex(WP), intent(in)              :: omega
    real(WP), intent(in)                 :: x(:)
    complex(WP), intent(in)              :: y(:,:)
    real(WP), intent(in)                 :: x_ref
    complex(WP), intent(in)              :: y_ref(:)
    type(ext_real_t), intent(in)         :: chi
    integer, intent(in)                  :: n_iter
    type(mode_t)                         :: md

    real(WP)    :: phase
    complex(WP) :: norm_fac
    integer     :: n_p
    integer     :: n_g
    integer     :: n_pg

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Construct the mode

    md%cf => cf

    md%op = op

    md%chi = chi

    md%x = x
    md%y = y

    md%x_ref = x_ref
    md%y_ref = y_ref

    md%omega = omega

    md%n = SIZE(md%x)
    md%n_iter = n_iter

    ! Normalize by the mode inertia, and so that y_ref(1) is real

    phase = ATAN2(AIMAG(md%y_ref(1)), REAL(md%y_ref(1)))

    norm_fac = 1._WP/SQRT(md%E())*EXP(CMPLX(0._WP, -phase, KIND=WP))

    md%y = norm_fac*md%y
    md%y_ref = norm_fac*md%y_ref

    ! Classify the mode

    call classify(md, n_p, n_g, n_pg)

    md%n_p = n_p
    md%n_g = n_g
    md%n_pg = n_pg

    ! Finish

    return

  end function init_md

!****

  $if($GFORTRAN_PR57922)

  subroutine final (this)

    class(mode_t), intent(inout) :: this

    ! Finalize the mode

    deallocate(this%x)
    deallocate(this%y)

    ! Finish

  end subroutine final

  $endif

!****

  $if($MPI)

  subroutine bcast_md (md, root_rank, cf)

    class(mode_t), intent(inout)        :: md
    integer, intent(in)                 :: root_rank
    class(coeffs_t), intent(in), target :: cf

    ! Broadcast the mode

    if(MPI_RANK /= root_rank) then
       md%cf => cf
    endif

    call bcast(md%op, root_rank)

    call bcast_alloc(md%x, root_rank)
    call bcast_alloc(md%y, root_rank)

    call bcast(md%omega, root_rank)

    call bcast(md%n_p, root_rank)
    call bcast(md%n_g, root_rank)
    call bcast(md%n_pg, root_rank)

    call bcast(md%n, root_rank)
    call bcast(md%n_iter, root_rank)

    ! Finish

    return

  end subroutine bcast_md

  $endif

!****

  function freq (this, freq_units)

    class(mode_t), intent(in)    :: this
    character(LEN=*), intent(in) :: freq_units
    complex(WP)                  :: freq

    ! Calculate the frequency

    freq = this%omega*freq_scale(this%cf, this%op, this%x(this%n), freq_units)

    ! Finish
    
    return

  end function freq

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

    associate (c_1 => this%cf%c_1(this%x), &
               l => this%op%l, omega_c => this%cf%omega_c(this%x, this%op%m, this%omega))

      if(l /= 0) then

         if(l /= 1) then

            where (this%x > 0._WP)
               xi_h = this%y(2,:)*this%x**(l-1)/(c_1*omega_c**2)
            elsewhere
               xi_h = 0._WP
            end where

         else

            xi_h = this%y(2,:)/(c_1*omega_c**2)

         endif

      else

         xi_h = 0._WP

      end if

    end associate

    ! Finish

    return

  end function xi_h

!****

  function Yt_1 (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: Yt_1(this%n)

    ! Calculate the Takata Y_1 function; this is based on [Tak2006,
    ! his eqn. 69]

    associate (J => 1._WP - this%cf%U(this%x)/3._WP)

      Yt_1 = J*this%y(1,:) + (this%y(3,:) - this%y(4,:))/3._WP
      
    end associate

    ! Finish

    return

  end function Yt_1

!****

  function Yt_2 (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: Yt_2(this%n)

    ! Calculate the Takata Y_2 function; this is based on [Tak2006,
    ! his eqn. 70], divided by V

    Yt_2 = this%y(2,:) - this%y(1,:) - this%y(3,:)

    ! Finish

  end function Yt_2

!****

  function phip (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: phip(this%n)
    
    ! Calculate the Eulerian gravitational potential perturbation in units of
    ! G M_star / R_star

    associate (c_1 => this%cf%c_1(this%x), l => this%op%l)

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

    associate (c_1 => this%cf%c_1(this%x), l => this%op%l)

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

  function delS_en (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delS_en(this%n)

    complex(WP) :: A_6(6,this%n)
    complex(WP) :: dy_6(this%n)
    complex(WP) :: y_5(this%n)

    ! Calculate the Lagrangian specific entropy perturbation in units
    ! of c_p, from the energy equation

    associate(V => this%cf%V(this%x), c_1 => this%cf%c_1(this%x), &
              nabla_ad => this%cf%nabla_ad(this%x), nabla => this%cf%nabla(this%x), &
              c_rad => this%cf%c_rad(this%x), dc_rad => this%cf%dc_rad(this%x), c_thm => this%cf%c_thm(this%x), &
              c_eps_ad => this%cf%c_eps_ad(this%x), c_eps_S => this%cf%c_eps_S(this%x), &              
              l => this%op%l, omega_c => this%cf%omega_c(this%x, this%op%m, this%omega))

      where(this%x /= 0)
         A_6(1,:) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
         A_6(2,:) = V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
         A_6(3,:) = l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad
         A_6(4,:) = 0._WP
         A_6(5,:) = c_eps_S - l*(l+1)*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
         A_6(6,:) = -1._WP - l
      elsewhere
         A_6(1,:) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
         A_6(2,:) = V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
         A_6(3,:) = l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad
         A_6(4,:) = 0._WP
         A_6(5,:) = -HUGE(0._WP)
         A_6(6,:) = -1._WP - l
      endwhere

      dy_6 = this%x*deriv(this%x, this%y(6,:))

      y_5 = (dy_6 - (A_6(1,:)*this%y(1,:) + &
                     A_6(2,:)*this%y(2,:) + &
                     A_6(3,:)*this%y(3,:) + &
                     A_6(4,:)*this%y(4,:) + &
                     A_6(6,:)*this%y(6,:)))/A_6(5,:)

      where(this%x /= 0._WP)
         delS_en = y_5*this%x**(l-2)
      elsewhere
         delS_en = 0._WP
      end where

    end associate
         
    ! Finish

    return

  contains

    function deriv (x, y) result (dy_dx)

      real(WP), intent(in)    :: x(:)
      complex(WP), intent(in) :: y(:)
      complex(WP)             :: dy_dx(SIZE(x))
      
      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Differentiate y(x) using centered finite differences

      n = SIZE(x)

      dy_dx(1) = (y(2) - y(1))/(x(2) - x(1))

      do i = 2,n-1
         dy_dx(i) = 0.5_WP*((y(i) - y(i-1))/(x(i) - x(i-1)) + &
                            (y(i+1) - y(i))/(x(i+1) - x(i)))
      end do

      dy_dx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))

      ! Finish

      return

    end function deriv

  end function delS_en

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

  function delL_rd (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delL_rd(this%n)

    complex(WP) :: A_5(6,this%n)
    complex(WP) :: dy_5(this%n)
    complex(WP) :: y_6(this%n)

    ! Calculate the Lagrangian luminosity perturbation in units of L_star, 
    ! from the radiative diffusion equation

    associate(V => this%cf%V(this%x), U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              nabla_ad => this%cf%nabla_ad(this%x), nabla => this%cf%nabla(this%x), &
              c_dif => this%cf%c_dif(this%x), c_rad => this%cf%c_rad(this%x), &
              kappa_S => this%cf%kappa_S(this%x), &
              l => this%op%l, omega_c => this%cf%omega_c(this%x, this%op%m, this%omega))

      A_5(1,:) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
      A_5(2,:) = V*(l*(l+1)/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
      A_5(3,:) = V*c_dif
      A_5(4,:) = V*nabla_ad
      A_5(5,:) = V*nabla*(4._WP - kappa_S) - (l - 2._WP)
      A_5(6,:) = -V*nabla/c_rad

      dy_5 = this%x*deriv(this%x, this%y(5,:))

      where(this%x /= 0._WP)

         y_6 = (dy_5 - (A_5(1,:)*this%y(1,:) + &
                        A_5(2,:)*this%y(2,:) + &
                        A_5(3,:)*this%y(3,:) + &
                        A_5(4,:)*this%y(4,:) + &
                        A_5(5,:)*this%y(5,:)))/A_5(6,:)

      elsewhere

         y_6 = 0._WP

      endwhere

      delL_rd = y_6*this%x**(l+1)

    end associate
         
    ! Finish

    return

  contains

    function deriv (x, y) result (dy_dx)

      real(WP), intent(in)    :: x(:)
      complex(WP), intent(in) :: y(:)
      complex(WP)             :: dy_dx(SIZE(x))
      
      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Differentiate y(x) using centered finite differences

      n = SIZE(x)

      dy_dx(1) = (y(2) - y(1))/(x(2) - x(1))

      do i = 2,n-1
         dy_dx(i) = 0.5_WP*((y(i) - y(i-1))/(x(i) - x(i-1)) + &
                            (y(i+1) - y(i))/(x(i+1) - x(i)))
      end do

      dy_dx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))

      ! Finish

      return

    end function deriv

  end function delL_rd

!****

  function delp (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: delp(this%n)

    ! Calculate the Lagrangian pressure perturbation in units of p

    associate (V => this%cf%V(this%x), pi_c => this%cf%pi_c(), l => this%op%l)

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

    associate (Gamma_1 => this%cf%Gamma_1(this%x))

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

    associate (nabla_ad => this%cf%nabla_ad(this%x))

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

    associate(U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l)
      dE_dx = (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*this%x**2/c_1
    end associate

    ! Finish

    return

  end function dE_dx

!****

  function dW_dx (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: dW_dx(this%n)

    ! Calculate the differential work in units of G M_star^2/R_star
    ! t_dyn/t_KH = t_dyn L_*.  The entropy-based expression for the
    ! work is used (cf. Unno et al.  1989, eqn. 25.9); the additional
    ! factor of 4 pi in the denominator comes from averaging over
    ! solid angle

    associate(c_thm => this%cf%c_thm(this%x))

      dW_dx = -PI*AIMAG(CONJG(this%delT())*this%delS())*c_thm*this%x**2/(4._WP*PI)

    end associate

    ! Finish

    return

  end function dW_dx

!****

  function I_0 (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: I_0(this%n)

    ! Calculate the I_0 integral (eqn. 42 of Takata 2006, PASJ, 58,
    ! 759). This should vanish for radial modes

    associate(U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l, x => this%x, y => this%y)

      I_0 = x**(l+1)*(U*y(1,:) + y(4,:))/c_1

    end associate

    ! Finish

    return

  end function I_0

!****

  function I_1 (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: I_1(this%n)

    ! Calculate the I_1 integral (eqn. 43 of Takata 2006, PASJ, 58,
    ! 759). This should vanish for dipole modes

    associate(U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l, omega_c => this%cf%omega_c(this%x, this%op%m, this%omega), &
              x => this%x, y => this%y)

      I_1 = x**(l+2)*(c_1*omega_c**2*U*y(1,:) - U*y(2,:) + &
                  (U - c_1*omega_c**2 - 2._WP)*y(3,:) + (c_1*omega_c**2 - 1._WP)*y(4,:))/c_1**2

    end associate

    ! Finish

    return

  end function I_1

!****

  function prop_type (this)

    class(mode_t), intent(in) :: this
    integer                   :: prop_type(this%n)

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    associate(x => this%x, V_g => this%cf%V(this%x)/this%cf%Gamma_1(this%x), &
              As => this%cf%As(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l, omega_c => REAL(this%cf%omega_c(this%x, this%op%m, this%omega)))

      prop_type = MERGE(1, 0, c_1*omega_c**2 > As) + &
                  MERGE(-1, 0, l*(l+1)/(c_1*omega_c**2) > V_g)

    end associate

    ! Finish

    return

  end function prop_type

!****

  function K (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: K(this%n)

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)

    ! Calculate the rotation splitting kernel

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(x => this%x, U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l)

      K = (ABS(xi_r)**2 + (l*(l+1)-1)*ABS(xi_h)**2 - 2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1

      K = K/integrate(x, K)

    end associate

    ! Finish

    return

  end function K

!****

  function beta (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: beta

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)

    ! Calculate the rotation splitting scale

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(x => this%x, U => this%cf%U(this%x), c_1 => this%cf%c_1(this%x), &
              l => this%op%l)

      beta = integrate(x, (ABS(xi_r)**2 + (l*(l+1)-1)*ABS(xi_h)**2 - 2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1) / &
             integrate(x, (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*x**2/c_1)

    end associate

    ! Finish

    return

  end function beta

!****

  function E (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E
    
    ! Calculate the total mode inertia [Aer2010, eqn. 3.139] in units
    ! of M_star R_star**2

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
    complex(WP) :: xi_r_1
    complex(WP) :: xi_h_1
    real(WP)    :: A2

    ! Calculate the normalized mode inertia (Aerts et al. 2010,
    ! eqn. 3.140)

    E = this%E()

    associate(l => this%op%l)

      select case (this%op%inertia_norm_type)
      case ('RADIAL')
         A2 = ABS(this%xi_r_ref())**2
      case ('HORIZ')
         A2 = l*(l+1)*ABS(this%xi_h_ref())**2
      case ('BOTH')
         A2 = ABS(this%xi_r_ref())**2 + l*(l+1)*ABS(this%xi_h_ref())**2
      case default
         $ABORT(Invalid inertia_norm_type)
      end select

      if(A2 == 0._WP) then
         $WARN(Amplitude at x_ref is zero; not normalizing inertia)
         E_norm = E
      else
         E_norm = E/A2
      endif

    end associate

    ! Finish

    return

  end function E_norm

!****

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

       associate(tau_thm => this%cf%tau_thm(this%x(i)), &
                 omega_c => this%cf%omega_c(this%x(i), this%op%m, this%omega))

         if(REAL(omega_c)*tau_thm/TWOPI > 1._WP) then
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

  function xi_r_ref (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: xi_r_ref
    
    ! Calculate the radial displacement perturbation at x_ref in units
    ! of R_star

    xi_r_ref = this%y_ref(1)

    ! Finish

    return

  end function xi_r_ref

!****

  function xi_h_ref (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: xi_h_ref
    
    ! Calculate the radial displacement perturbation at x_ref in units
    ! of R_star

    associate (c_1 => this%cf%c_1(this%x_ref), &
               l => this%op%l, omega_c => this%cf%omega_c(this%x_ref, this%op%m, this%omega))

      if(l /= 0) then

         xi_h_ref = this%y_ref(2)/(c_1*omega_c**2)

      else

         xi_h_ref = 0._WP

      end if

    end associate

    ! Finish

    return

  end function xi_h_ref

!****

  subroutine classify (md, n_p, n_g, n_pg)

    class(mode_t), intent(in) :: md
    integer, intent(out)      :: n_p
    integer, intent(out)      :: n_g
    integer, intent(out)      :: n_pg

    real(WP) :: y_1(md%n)
    real(WP) :: y_2(md%n)
    real(WP) :: x_turn
    integer  :: i
    integer  :: n_c
    integer  :: n_a

    ! Classify the mode based on its eigenfunctions

    select case (md%op%l)

    case (0)

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       y_1 = REAL(md%y(1,:))
       y_2 = REAL(md%y(2,:))

       mono_loop : do i = 2,md%n-1
          if((y_1(i) >= y_1(i-1) .AND. y_1(i+1) >= y_1(i)) .OR. &
             (y_1(i) <= y_1(i-1) .AND. y_1(i+1) <= y_1(i))) exit mono_loop
       end do mono_loop

       ! Count winding numbers

       call count_windings(y_1(i:), y_2(i:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       n_p = n_a + 1
       n_g = n_c

       n_pg = n_p - n_g

    case (1)

       ! Dipole modes

       ! Set up the Takata Y^a_1 and Y^a_2 functions
       
       y_1 = REAL(md%Yt_1())
       y_2 = REAL(md%Yt_2())

       ! Find the inner turning point (this is to deal with noisy
       ! near-zero solutions at the origin)

       call find_x_turn(md%x, md%cf, md%op, REAL(md%omega), x_turn)

       x_turn_loop : do i = 1,md%n-1
          if(md%x(i) > x_turn) exit x_turn_loop
       end do x_turn_loop

       ! Count winding numbers

!       call count_windings(y_1(i:), y_2(i:), n_c, n_a, md%x)
       call count_windings(y_1(i:), y_2(i:), n_c, n_a)

       n_p = n_a
       n_g = n_c

       if(n_p >= n_g) then
          n_pg = n_p - n_g + 1
       else
          n_pg = n_p - n_g
       endif

    case default

       ! Other modes

       y_1 = REAL(md%y(1,:))
       y_2 = REAL(md%y(2,:))

       ! Count winding numbers

       call count_windings(y_1, y_2, n_c, n_a)

       ! Classify

       n_p = n_a
       n_g = n_c

       n_pg = n_p - n_g

    end select

    ! Finish

    return

  contains

    subroutine count_windings (y_1, y_2, n_c, n_a, x)

      real(WP), intent(in)           :: y_1(:)
      real(WP), intent(in)           :: y_2(:)
      integer, intent(out)           :: n_c
      integer, intent(out)           :: n_a
      real(WP), intent(in), optional :: x(:)

      integer  :: i
      real(WP) :: y_2_cross

      $CHECK_BOUNDS(SIZE(y_2),SIZE(y_1))

      if(PRESENT(x)) then
         $CHECK_BOUNDS(SIZE(x),SIZE(y_1))
      endif

      ! Count clockwise (n_c) and anticlockwise (n_a) windings in the (y_1,y_2) plane

      n_c = 0
      n_a = 0

      do i = 1,SIZE(y_1)-1

         ! Look for a node in y_1

         if(y_1(i) >= 0._WP .AND. y_1(i+1) < 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

            if(y_2_cross >= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(i),x(i+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(i),x(i+1)
            endif

         elseif(y_1(i) <= 0._WP .AND. y_1(i+1) > 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

            if(y_2_cross <= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(i),x(i+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(i),x(i+1)
            endif

         endif

      end do

      ! Finish

      return

    end subroutine count_windings

  end subroutine classify

end module gyre_mode
