! Module   : gyre_sol
! Purpose  : solution data
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
$include 'core_memory.inc'

module gyre_sol

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants
  use gyre_ext
  use gyre_freq
  use gyre_grid
  use gyre_model
  $if ($MPI)
  use gyre_model_mpi
  $endif
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure       :: ${NAME}_1_
    procedure       :: ${NAME}_n_
    generic, public :: ${NAME} => ${NAME}_1_, ${NAME}_n_
  $endsub

  type :: sol_t
     private
     class(model_t), pointer, public  :: ml => null()
     class(c_rot_t), allocatable      :: rt
     real(WP), allocatable, public    :: x(:)
     integer, allocatable, public     :: s(:)
     complex(WP), allocatable, public :: y(:,:)
     complex(WP)                      :: omega
     complex(WP)                      :: l_i
     type(c_ext_t)                    :: discrim
     type(r_ext_t), public            :: chi
     integer                          :: k_ref
     integer, public                  :: n_pg
     integer, public                  :: n_p
     integer, public                  :: n_g
     integer, public                  :: n_iter
     integer, public                  :: n_k
     logical                          :: pruned
   contains
     private
     procedure, public :: prune => prune_
     procedure, public :: freq => freq_
     $PROC_DECL(xi_r)
     $PROC_DECL(xi_h)
     $PROC_DECL(eul_phi)
     $PROC_DECL(deul_phi)
     $PROC_DECL(lag_S)
     $PROC_DECL(lag_L)
     $PROC_DECL(eul_P)
     $PROC_DECL(lag_P)
     $PROC_DECL(eul_rho)
     $PROC_DECL(lag_rho)
     $PROC_DECL(eul_T)
     $PROC_DECL(lag_T)
     $PROC_DECL(dE_dx)
     $PROC_DECL(dW_dx)
     $PROC_DECL(F_j)
     $PROC_DECL(Yt_1)
     $PROC_DECL(Yt_2)
     $PROC_DECL(I_0)
     $PROC_DECL(I_1)
     procedure, public :: lag_T_eff => lag_T_eff_
     procedure, public :: lag_g_eff => lag_g_eff_
     procedure, public :: f_T => f_T_
     procedure, public :: f_g => f_g_
     procedure, public :: psi_T => psi_T_
     procedure, public :: E => E_
     procedure, public :: E_p => E_p_
     procedure, public :: E_g => E_g_
     procedure, public :: E_norm => E_norm_
     procedure, public :: W => W_
     procedure, public :: K => K_
     procedure, public :: beta => beta_
     procedure, public :: omega_int => omega_int_
     procedure, public :: eta => eta_
     procedure, public :: prop_type => prop_type_
     procedure         :: classify_
  end type mode_t

  ! Interfaces

  interface sol_t
     module procedure sol_t_
  end interface sol_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: sol_t
  public :: reallocate
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function sol_t_ (ml, s, x, y, omega, k_ref, md_p, os_p) result (sl)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s(:)
    real(WP), intent(in)                :: x(:)
    complex(WP), intent(in)             :: y(:,:)
    complex(WP), intent(in)             :: omega
    integer, intent(in)                 :: k_ref
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(sol_t)                         :: sl

    integer  :: k
    real(WP) :: phase

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(s))

    ! Construct the sol_t

    sl%ml => ml

    allocate(sl%co, SOURCE=co)
    allocate(sl%y, SOURCE=y)

    sl%n_k = SIZE(co)

    allocate(sl%rt, SOURCE=c_rot_t(ms, md_p, os_p))

    sl%omega = omega
    sl%l_i = sl%rt%l_i(omega)

    sl%k_ref = k_ref
    
    ! Normalize the solution so that y_ref(1) is purely real, and the
    ! total mode energy is unity

    phase = ATAN2(AIMAG(sl%y(1,sl%k_ref)), &
                  REAL(sl%y(1,sl%k_ref)))

    sl%y = sl%y/SQRT(sl%E())*EXP(CMPLX(0._WP, -phase, KIND=WP))

    ! Classify the solution

    call this%classify_()

    ! Do other initializations

    sl%md_p = md_p
    sl%os_p = os_p

    sl%pruned = .FALSE.

    ! Finish

    return

  end function sol_t_

  !****

  $REALLOCATE(type(sol_t),1)

  !****

  subroutine prune_ (this)

    class(sol_t), intent(inout) :: this

    ! Prune the sol_t

    if (.NOT. this%pruned) then

       deallocate(this%s)
       deallocate(this%x)
       deallocate(this%y)

       this%n_k = 0

       this%pruned = .TRUE.

    endif

    ! Finish

    return

  end subroutine prune_

  !****

  function freq_ (this, freq_units, freq_frame) result (freq)

    class(sol_t), intent(in)           :: this
    character(*), intent(in)           :: freq_units
    character(*), optional, intent(in) :: freq_frame
    complex(WP)                        :: freq

    ! Calculate the frequency

    if (PRESENT(freq_frame)) then
       freq = freq_from_omega(this%omega, this%ml, this%mp, this%op, this%x(1), this%x(this%n), freq_units, freq_frame)
    else
       freq = freq_from_omega(this%omega, this%ml, this%mp, this%op, this%x(1), this%x(this%n), freq_units, 'INERTIAL')
    endif

    ! Finish
    
    return

  end function freq_

  !****

  function xi_r_1_ (this, k) result (xi_r)

    class(sol_seg_t), intent(in) :: this
    integer, intent(in)          :: k
    complex(WP)                  :: xi_r

    ! Evaluate the radial displacement perturbation, in units of
    ! R_star

    associate (x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      if (l_i /= 1._WP) then

         if (x /= 0._WP) then
            xi_r = y(1)*x**(l_i-1._WP)
         else
            xi_r = 0._WP
         endif

      else

         xi_r = y(1)

      endif

    end associate

    ! Finish

    return

  end function xi_r_1_

  !****

  function xi_h_1_ (this, k) result (xi_h)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: xi_h

    real(WP)    :: c_1
    complex(WP) :: omega_c

    ! Evaluate the horizontal displacement perturbation, in units of
    ! R_star

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      c_1 = this%ml%c_1(x, s)

      omega_c = this%rt%omega_c(x, s, this%omega)
      
      if (l_i /= 0._WP) then

         if (l_i /= 1._WP) then

            if (x /= 0._WP) then
               xi_h = y(2)*x**(l_i-1._WP)/(c_1*omega_c**2)
            else
               xi_h = 0._WP
            end if

         else
            
            xi_h = y(2)/(c_1*omega_c**2)

         endif

      else

         xi_h = 0._WP

      end if

    end associate

    ! Finish

    return

  end function xi_h_1_

  !****

  function eul_phi_1_ (this, k) result (eul_phi)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: eul_phi

    real(WP) :: c_1
    
    ! Evaluate the Eulerian gravitational potential perturbation, in
    ! units of G M_star / R_star

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i = this%l_i)

      c_1 = this%ml%c_1(s, x)

      if (l_i /= 0._WP) then

         if (x /= 0._WP) then
            eul_phi = y(3)*x**l_i/c_1
         else
            eul_phi = 0._WP
         endif

      else

         eul_phi = y(3)/c_1

      endif

    end associate

    ! Finish

    return

  end function eul_phi_1_

  !****

  function deul_phi_1_ (this, k) result (deul_phi)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: deul_phi

    real(WP) :: c_1
    
    ! Evaluate the Eulerian potential gradient (gravity) perturbation,
    ! in units of G M_star / R_star**2

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      c_1 = this%ml%c_1(s, x)

      if (l_i /= 1._WP) then

         if (x /= 0._WP) then
            deul_phi = y(4)*x**(l_i-1._WP)/c_1
         else
            deul_phi = 0._WP
         end if
       
      else

         deul_phi = y(4)/c_1

      end if

    end associate

    ! Finish

    return

  end function deul_phi_1_

  !****

  function lag_S_1_ (this, k) result (lag_S)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lag_S

    ! Evaluate the Lagrangian specific entropy perturbation, in units
    ! of c_p

    associate (x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      if (x /= 0._WP) then
         lag_S = y(5)*x**(l_i-2._WP)
      else
         lag_S = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_S_1_

  !****

  function lag_L_1_ (this, k) result (lag_L)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lag_L

    ! Evaluate the Lagrangian radiative luminosity perturbation, in
    ! units of L_star

    associate (x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      if (x /= 0._WP) then
         lag_L = y(6)*x**(l_i+1._WP)
      else
         lag_L = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_L_1_

  !****

  function eul_P_1_ (this, k) result (eul_P)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: eul_P

    complex(WP) :: xi_r
    complex(WP) :: lag_P
    real(WP)    :: V_2

    ! Evaluate the Eulerian pressure perturbation, in units of P

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      lag_P = this%lag_P(k)

      V_2 = this%ml%V_2(s, x)

      eul_P = lag_P + V_2*x*xi_r

    end associate

    ! Finish

    return

  end function eul_P_1_
  
  !****

  function lag_P_1_ (this, k) result (lag_P)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lag_P

    real(WP) :: V_2

    ! Evaluate the Lagrangian pressure perturbation, in units of P

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      V_2 = this%ml%V_2(s, x)

      if (l_i /= 0._WP) then

         if (x /= 0._WP) then
            lag_P = V_2*(y(2) - y(1) - y(3))*x**l_i
         else
            lag_P = 0._WP
         end if

      else

         lag_P = V_2*(y(2) - y(1) - y(3))

      endif

    end associate

    ! Finish

    return

  end function lag_P_1_

  !****

  function eul_rho_1_ (this, k) result (eul_rho)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: eul_rho

    complex(WP) :: xi_r
    complex(WP) :: lag_rho
    real(WP)    :: U
    real(WP)    :: dU
    real(WP)    :: D

    ! Evaluate the Eulerian density perturbation, in units of rho

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      lag_rho = this%lag_rho(k)

      U = this%ml%U(s, x)
      dU = this%ml%dU(s, x)

      D = dU + U - 3._WP

      if (x /= 0._WP) then
         eul_rho = lag_rho - D*xi_r/x
      else
         eul_rho = lag_rho
      endif

    end associate

    ! Finish

    return

  end function eul_rho_1_

  !****

  function lag_rho_1_ (this, k) result (lag_rho)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lag_rho

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: Gamma_1
    real(WP)    :: delta

    ! Evaluate the Lagrangian density perturbation, in units of
    ! rho. This expression implements eqn. 13.83 of [Unn1989]

    associate (s => this%s(k), &
               x => this%x(k))

      lag_P = this%lag_P(k)
      lag_S = this%lag_S(k)

      Gamma_1 = this%ml%Gamma_1(s, x)
      delta = this%ml%delta(s, x)

      lag_rho = lag_P/Gamma_1 - delta*lag_S

    end associate

    ! Finish

    return

  end function lag_rho_1_

  !****

  function eul_T_1_ (this, k) result (eul_T)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: eul_T

    complex(WP) :: xi_r
    complex(WP) :: lag_T
    real(WP)    :: V_2
    real(WP)    :: nabla

    ! Evaluate the Lagrangian temperature perturbation, in units of T

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      lag_T = this%lag_T(k)

      V_2 = this%ml%V_2(s, x)
      nabla = this%ml%nabla(s, x)
      
      eul_T = lag_T + nabla*V_2*x*xi_r

    end associate

    ! Finish

    return

  end function eul_T_1_

  !****

  function lag_T_1_ (this, k) result (lag_T)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lag_T

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: nabla_ad

    ! Evaluate the Lagrangian temperature perturbation, in units of
    ! T. This expression implements eqn. 13.84 of [Unn1989]

    associate (s => this%s(k), &
               x => this%x(k))

      lag_P = this%lag_P(k)
      lag_S = this%lag_S(k)

      nabla_ad = this%ml%nabla_ad(s, x)
      
      lag_T = nabla_ad*lag_P + lag_S

    end associate

    ! Finish

    return

  end function lag_T_1_

  !****

  function lambda_1_ (this, k) result (lambda)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: lambda

    ! Evaluate the angular eigenvalue

    associate (this%s => this%s(k), &
               this%x => this%x(k))

      lambda = this%rt%lambda(s, x, this%omega)

    end associate

    ! Finish

    return

  end function lambda_1_
    
  !****

  function dE_dx_1_ (this, k) result (dE_dx)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    real(WP)                 :: dE_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: U
    real(WP)    :: c_1

    ! Evaluate the differential mode inertia, in units of M_star
    ! R_star**2. This expression is based on eqn. 3.139 of [Aer2010],
    ! with the initial factor of 4 pi cancelled by their definitions of
    ! \tilde{\xi}_r and \tilde{\xi}_h (cf. eqns. 3.124 and 3.131, ibid)

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      lambda = this%lambda(k)

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      dE_dx = (ABS(xi_r)**2 + ABS(lambda)*ABS(xi_h)**2)*U*x**2/(4._WP*PI*c_1)

    end associate

    ! Finish

    return

  end function dE_dx_1_

  !****

  function dW_dx_1_ (this, k) result (dW_dx)

    use gyre_evol_model

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    real(WP)                 :: dW_dx

    real(WP)    :: t_dyn
    real(WP)    :: t_kh
    complex(WP) :: lag_T
    complex(WP) :: lag_S
    real(WP)    :: c_thm

    ! Evaluate the differential work, in units of G M_star**2/R_star.
    ! This expression is based on eqn. 25.9 of [Unn1989]

    select type (ml => this%ml)
    class is (evol_model_t)
       t_dyn = SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    associate (s => this%s(k), &
               x => this%x(k))

      lag_T = this%lag_T(k)
      lag_S = this%lag_S(k)
    
      c_thm = this%ml%c_thm(s, x)

      dW_dx = PI*AIMAG(CONJG(lag_T)*lag_S)*c_thm*x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_dx_1_

  !****

  function K_n_1_ (this, k) result (K_n)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    real(WP)                 :: K_n

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    integer     :: l
    real(WP)    :: U
    real(WP)    :: c_1
    
    ! Calculate the numerator of the rotation splitting kernel. This
    ! expression is based on equation 3.356 of [Aer2010]

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      l = this%rt%l

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      K_n = (ABS(xi_r)**2 + (l*(l+1)-1)*ABS(xi_h)**2 - &
             2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1

    end associate

    ! Finish

    return

  end function K_n_1_

  !****

  function F_j_1_ (this, k) result (F_j)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    real(WP)                 :: F_j

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: c_1
    real(WP)    :: U
    integer     :: m
    complex(WP) :: omega_c
    
    ! Evaluate the angle-averaged angular momentum flux due to
    ! Reynolds stress, in units of G M_star**2/R_star**3.  This
    ! expression is based on eqn. 21 of [LeeSai1993]

    associate (s => this%s(k), &
               x => this%x(k))

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      c_1 = this%ml%c_1(s, x)
      U = this%ml%U(s, x)

      m = this%rt%m

      omega_c = this%rt%omega_c(s, x, this%omega)

      F_j = -m*ABS(omega_c**2)*x*U*AIMAG(CONJG(xi_r)*xi_h)/(32._WP*PI**2*c_1)

    end associate

    ! Finish

    return

  end function F_j_1_

  !****

  function Yt_1_1_ (this, k) result (Yt_1)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: Yt_1

    real(WP) :: J

    ! Evaluate the Takata Y_1 function. This expression is based on
    ! eqn. 69 of [Tak2006b]

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k))

      J = 1._WP - this%ml%U(s, x)/3._WP

      Yt_1 = J*y(1) + (y(3) - y(4))/3._WP

    end associate
      
    ! Finish

    return

  end function Yt_1_1_

  !****

  function Yt_2_1_ (this, k) result (Yt_2)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: Yt_2

    ! Evaluate the Takata Y_2 function. This expression is based on
    ! eqn. 70 of [Tak2006b], divided through by V

    associate (y => this%y(:,k))

      Yt_2 = y(2) - y(1) - y(3)

    end associate

    ! Finish

    return

  end function Yt_2_1_

  !****

  function I_0_1_ (this, k) result (I_0)

    class(sol_seg_t), intent(in) :: this
    integer, intent(in)          :: k
    complex(WP)                  :: I_0

    real(WP) :: U
    real(WP) :: c_1
    
    ! Evaluate the I_0 integral, which should be zero for radial
    ! modes. This expression is based on eqn. 42 of [Tak2006a]

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      U = this%ml%U(s, x)
      c_1 = this%ms%c_1(s, x)

      if (x /= 0._WP) then
         I_0 = x**(l_i+1._WP)*(U*y(1) + y(4))/c_1
      else
         I_0 = 0._WP
      endif

    end associate

    ! Finish

    return

  end function I_0_1_

  !****

  function I_1_1_ (this, k) result (I_1)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: k
    complex(WP)              :: I_1

    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: omega_c

    ! Evaluate the I_0 integral, which should be zero for dipole
    ! modes. This expression is based on eqn. 43 of [Tak2006a]

    associate (s => this%s(k), &
               x => this%x(k), &
               y => this%y(:,k), &
               l_i => this%l_i)

      U = this%ml%U(s, x)
      c_1 = this%ml%U(s, x)

      omega_c = this%rt%omega_c(s, x, this%omega)

      if (x /= 0._WP) then
         I_1 = x**(l_i+2._WP)*(c_1*omega_c**2*U*y(1) - U*y(2) + &
               (U - c_1*omega_c**2 - 2._WP)*y(3) + (c_1*omega_c**2 - 1._WP)*y(4))/c_1**2
      else
         I_1 = 0._WP
      endif

    end associate
      
    ! Finish

    return

  end function I_1_1_

  !****

  $define $PROC_N $sub

  $local $NAME $1

  function ${NAME}_n_ (this) result (${NAME})

    class(sol_t), intent(in) :: this
    complex(WP)              :: $NAME(this%n_k)

    integer :: k

    ! Evaluate $NAME

    !$OMP PARALLEL DO
    k_loop : do k = 1, this%n
       ${NAME}(k) = this%${NAME}(k)
    end do k_loop

    ! Finish

    return

  end function ${NAME}_n_

  $endsub

  $PROC_N(xi_r)
  $PROC_N(xi_h)
  $PROC_N(eul_phi)
  $PROC_N(deul_phi)
  $PROC_N(lag_S)
  $PROC_N(lag_L)
  $PROC_N(eul_P)
  $PROC_N(lag_P)
  $PROC_N(eul_rho)
  $PROC_N(lag_rho)
  $PROC_N(eul_T)
  $PROC_N(lag_T)
  $PROC_N(dE_dx)
  $PROC_N(dW_dx)
  $PROC_N(F_j)
  $PROC_N(Yt_1)
  $PROC_N(Yt_2)
  $PROC_N(I_0)
  $PROC_N(I_1)

  !****

  function lag_T_eff_ (this) result (lag_T_eff)

    class(sol_t), intent(in) :: this
    complex(WP)              :: lag_T_eff

    complex(WP) :: xi_r
    complex(WP) :: lag_L

    ! Calculate the effective temperature perturbation at x_ref
    ! (assumed to correspond to the photosphere), in units of
    ! T_eff. This expression is based on the standard definition of
    ! effective temperature

    xi_r = this%xi_r(this%k_ref)
    lag_L = this%lag_L(this%k_ref)

    lag_T_eff = 0.25_WP*(lag_L - 2._WP*xi_r)

    ! Finish

    return

  end function lag_T_eff_

  !****

  function lag_g_eff_ (this) result (lag_g_eff)

    class(sol_t), intent(in) :: this
    complex(WP)              :: lag_g_eff

    complex(WP) :: xi_r
    complex(WP) :: deul_phi
    real(WP)    :: c_1
    real(WP)    :: U

    ! Calculate the effective gravity perturbation at x_ref (assumed
    ! to correspond to the photosphere), in units of the gravity. This
    ! expression is based on eqn. 24 of [Dup2002]

    associate (s => this%s(this%k_ref), &
               x => this%x(this%k_ref), &
               omega => this%omega)

      xi_r = this%xi_r(this%k_ref))
      deul_phi = this%deul_phi(this%k_ref)

      c_1 = this%ml%c_1(s, x)
      U = this%ml%c_1(s, x)

      lag_g_eff = (c_1/x)*deul_phi + (U - (2._WP + c_1*omega**2))*xi_r/x

    end associate

    ! Finish

    return

  end function lag_g_eff_

  !****

  function f_T_ (this) result (f_T)

    class(sol_t), intent(in) :: this
    real(WP)                 :: f_T

    complex(WP) :: C_T

    ! Evaluate the non-adiabatic f_T parameter. This is expression is
    ! based on eqn. 5 of [Dup2003]

    C_T = this%lag_T_eff()/this%xi_r(this%k_ref)

    f_T = ABS(C_T)

    ! Finish

    return

  end function f_T_

  !****

  function f_g_ (this) result (f_g)

    class(sol_t), intent(in) :: this
    real(WP)                 :: f_g

    complex(WP) :: C_g

    ! Evaluate the non-adiabatic f_g parameter. This is expression is
    ! based on eqn. 6 of [Dup2003]

    C_g = this%lag_g_eff()/this%xi_r(this%k_ref)

    f_g = -ABS(C_g)

    ! Finish

    return

  end function f_g_

  !****

  function psi_T_ (this) result (psi_T)

    class(sol_t), intent(in) :: this
    real(WP)                 :: psi_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic psi_T parameter, in radians. This is
    ! expression is based on eqn. 5 of [Dup2003]

    C_T = this%lag_T_eff()/this%xi_r(this%k_ref)

    psi_T = ATAN2(AIMAG(C_T), REAL(C_T))

    ! Finish

    return

  end function psi_T_

  !****

  function E_ (this) result (E)

    class(sol_t), intent(in) :: this
    real(WP)                 :: E
    
    ! Calculate the mode inertia, in units of M_star R_star**2

    E = integrate(this%x, this%dE_dx())

    ! Finish

    return

  end function E_

  !****

  function E_p_ (this) result (E_p)

    class(sol_t), intent(in) :: this
    real(WP)                 :: E_p

    ! Calculate the mode inertia in acoustic-wave propagation regions,
    ! in units of M_star R_star**2

    E_p = integrate(this%x, this%dE_dx(), mask=(this%prop_type == 1))
    
    ! Finish

    return

  end function E_p_

  !****

  function E_g_ (this) result (E_g)

    class(sol_t), intent(in) :: this
    real(WP)                 :: E_g

    ! Calculate the mode inertia in gravity-wave propagation regions,
    ! in units of M_star R_star**2

    E_g = integrate(this%x, this%dE_dx(), mask=(this%prop_type == -1))
    
    ! Finish

    return

  end function E_g_

  !****

  function E_norm_ (this) result (E_norm)
 
    class(sol_t), intent(in) :: this
    real(WP)                 :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: A2

    ! Calculate the normalized mode inertia. This expression is based
    ! on eqn. 3.140 of [Aer2010]

    E = this%E()

    xi_r = this%xi_r(this%k_ref)
    xi_h = this%xi_h(this%k_ref))

    lambda = this%lambda(this%k_ref)

    select case (this%op%inertia_norm)
    case ('RADIAL')
       A2 = ABS(xi_r)**2
    case ('HORIZ')
       A2 = ABS(lambda)*ABS(xi_h)**2
    case ('BOTH')
       A2 = ABS(xi_r)**2 + ABS(lambda)*ABS(xi_h)**2
    case default
       $ABORT(Invalid inertia_norm)
    end select

    if (A2 == 0._WP) then
       $WARN(Amplitude at x_ref is zero; not normalizing inertia)
       E_norm = E
    else
       E_norm = E/A2
    endif

    ! Finish

    return

  end function E_norm_

  !****

  function W_ (this) result (W)

    class(sol_t), intent(in) :: this
    real(WP)                 :: W
    
    ! Calculate the total work, in units of G M_star**2/R_star

    W = integrate(this%x, this%dW_dx())

    ! Finish

    return

  end function W_

  !****

  function K_ (this) result (K)

    class(sol_t), intent(in) :: this
    complex(WP)              :: K(this%n_k)

    real(WP) :: K_n(this%n_k)
    real(WP) :: K_d

    ! Calculate the rotation splitting kernel

    K_n = this%K_n()
    K_d = integrate(this%x, K_n)

    K = K_n/K_d

    ! Finish

    return

  end function K_

  !****

  function beta_ (this) result (beta)

    class(sol_t), intent(in) :: this
    complex(WP)              :: beta

    complex(WP) :: xi_r(this%n_k)
    complex(WP) :: xi_h(this%n_k)
    integer     :: l
    real(WP)    :: U(this%n_k)
    real(WP)    :: c_1(this%c_1)
    real(WP)    :: dE_dx_0(this%n_k)

    ! Calculate the rotation splitting scale

    associate (s => this%s, &
               x => this%x)

      xi_r = this%xi_r()
      xi_h = this%xi_h()

      l = this%rt%l

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      dE_dx_0 = (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*x**2/(4._WP*PI*c_1)

    end associate

    beta = integrate(this%x, this%K_d())/ &
           integrate(this%x, dE_dx_0)
      
    ! Finish

    return

  end function beta_

  !****

  function omega_int_ (this) result (omega_int)

    class(sol_t), intent(in) :: this
    complex(WP)              :: omega_int

    complex(WP) :: xi_r(this%n_k)
    complex(WP) :: eul_phi(this%n_k)
    complex(WP) :: eul_rho(this%n_k)
    complex(WP) :: lag_rho(this%n_k)
    complex(WP) :: lag_P(this%n_k)
    real(WP)    :: V_2(this%n_k)
    real(WP)    :: As(this%n_k)
    real(WP)    :: U(this%n_k)
    real(WP)    :: c_1(this%n_k)
    real(WP)    :: Gamma_1(this%n_k)
    real(WP)    :: V_g(this%n_k)
    real(WP)    :: x4_V(this%n_k)
    complex(WP) :: W_th
    complex(WP) :: W_re
    complex(WP) :: W_gr
    complex(WP) :: W_xi

    ! Calculate the dimensionless frequency from the integral
    ! expression in eqn. (1.71) of [Dup2003]

    associate (s => this%s, &
         x => this%x)

      xi_r = this%xi_r()
      eul_phi = this%eul_phi()
      eul_rho = this%eul_rho()
      lag_rho = this%lag_rho()
      lag_P = this%lag_P()

      V_2 = this%ml%V_2(s, x)
      As = this%ml%As(s, x)
      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      Gamma_1 = this%ml%Gamma_1(s, x)

      V_g = V_2*x**2/Gamma_1
      x4_V = x**2/V_2

      W_th = integrate(x, CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2)))

      W_re = integrate(x, 2._WP*REAL(lag_rho*CONJG(xi_r)*(x/c_1)*(x**2*U/c_1)))

      W_gr = integrate(x, CONJG(eul_rho)*eul_phi*(x**2*U/c_1))

      W_xi = integrate(x, -ABS(xi_r)**2*(x/c_1)*(x*U*(-V_g-As)/c_1))

      omega_int = SQRT(4._WP*PI*(W_th + W_re + W_gr + W_xi)/this%E())

    end associate
    
    ! Finish

    return

  end function omega_int_

  !****

  function eta_ (this) result (eta)

    class(mode_t), intent(in) :: this
    real(WP)                  :: eta

    real(WP) :: dW_dx(this%n_k)

    ! Calculate the normalized growth rate defined by [Stel1978]

    dW_dx = this%dW_dx()

    eta = integrate(this%x, dW_dx)/integrate(this%x, ABS(dW_dx))

    ! Finish

    return

  end function eta_
    
  !****

  function prop_type_ (this) result (prop_type)

    class(mode_t), intent(in) :: this
    integer                   :: prop_type(this%n_k)

    real(WP)    :: V_2(this%n_k)
    real(WP)    :: As(this%n_k)
    real(WP)    :: c_1(this%n_k)
    real(WP)    :: Gamma_1(this%n_k)
    complex(WP) :: lambda(this%n_k)
    real(WP)    :: V_g(this%n_k)

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    associate (s => this%s, &
               x => this%x)

      V_2 = this%ml%V_2(s, x)
      As = this%ml%As(s, x)
      c_1 = this%ml%c_1(s, x)

      Gamma_1 = this%ml%Gamma_1(s, x)

      lambda = this%lambda()

      V_g = V_2*x**2/Gamma_1

      prop_type = MERGE(1, 0, REAL(c_1*omega_c**2) > As) + &
                  MERGE(-1, 0, REAL(lambda/(c_1*omega_c**2)) > V_g)

    end associate

    ! Finish

    return

  end function prop_type_

  !****

  subroutine classify_ (this)

    class(sol_t), intent(inout) :: this

    real(WP) :: y_1(this%n_k)
    real(WP) :: y_2(this%n_k)
    real(WP) :: x_turn
    integer  :: k
    integer  :: n_c
    integer  :: n_a

    ! Classify the mode based on its eigenfunctions

    l = this%rt%l

    if (l == 0) then

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       y_1 = REAL(this%y(1,:))
       y_2 = REAL(this%y(2,:))

       mono_loop : do k = 2, this%n_k-1
          if ((y_1(k) >= y_1(k-1) .AND. y_1(k+1) >= y_1(k)) .OR. &
              (y_1(k) <= y_1(k-1) .AND. y_1(k+1) <= y_1(k))) exit mono_loop
       end do mono_loop

       ! Count winding numbers

       call count_windings_(y_1(k:), y_2(k:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       this%n_p = n_a + 1
       this%n_g = n_c

       this%n_pg = this%n_p - this%n_g

    elseif (l == 1 .AND. .NOT. this%os_p%cowling_approx) then

       ! Dipole modes (non-Cowling)

       ! Set up the Takata Y^a_1 and Y^a_2 functions
       
       y_1 = REAL(this%Yt_1())
       y_2 = REAL(this%Yt_2())

       ! Find the inner turning point (this is to deal with noisy
       ! near-zero solutions at the origin)

       call find_x_turn(this%x, this%s, this%ml, this%md_p, this%os_p, REAL(this%omega), x_turn)

       x_turn_loop : do k = 1, this%n_k-1
          if (this%x(k) > x_turn) exit x_turn_loop
       end do x_turn_loop

       ! Count winding numbers

!       call count_windings(y_1(i:), y_2(i:), n_c, n_a, this%x)
       call count_windings_(y_1(i:), y_2(i:), n_c, n_a)

       this%n_p = n_a
       this%n_g = n_c

       if (this%n_p >= this%n_g) then
          this%n_pg = this%n_p - this%n_g + 1
       else
          this%n_pg = this%n_p - this%n_g
       endif

    else

       ! Other modes

       y_1 = REAL(this%y(1,:))
       y_2 = REAL(this%y(2,:))

       ! Count winding numbers

       call count_windings_(y_1, y_2, n_c, n_a)

       ! Classify

       this%n_p = n_a
       this%n_g = n_c

       this%n_pg = this%n_p - this%n_g

    endif

    ! Finish

    return

  contains

    subroutine count_windings_ (y_1, y_2, n_c, n_a, x)

      real(WP), intent(in)           :: y_1(:)
      real(WP), intent(in)           :: y_2(:)
      integer, intent(out)           :: n_c
      integer, intent(out)           :: n_a
      real(WP), optional, intent(in) :: x(:)

      integer  :: k
      real(WP) :: y_2_cross

      $CHECK_BOUNDS(SIZE(y_2),SIZE(y_1))

      if(PRESENT(x)) then
         $CHECK_BOUNDS(SIZE(x),SIZE(y_1))
      endif

      ! Count clockwise (n_c) and anticlockwise (n_a) windings in the (y_1,y_2) plane

      n_c = 0
      n_a = 0

      do k = 1,SIZE(y_1)-1

         ! Look for a node in y_1

         if (y_1(k) >= 0._WP .AND. y_1(k+1) < 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(k) - y_1(k)*(y_2(k+1) - y_2(k))/(y_1(k+1) - y_1(k))

            if(y_2_cross >= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(k),x(k+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(k),x(k+1)
            endif

         elseif (y_1(k) <= 0._WP .AND. y_1(k+1) > 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(k) - y_1(k)*(y_2(k+1) - y_2(k))/(y_1(k+1) - y_1(k))

            if (y_2_cross <= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(k),x(k+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(k),x(k+1)
            endif

         endif

      end do

      ! Finish

      return

    end subroutine count_windings_

  end subroutine classify_

  !****

  $if($MPI)

  subroutine bcast_ (sl, root_rank, ml)

    class(mode_t), intent(inout)       :: md
    integer, intent(in)                :: root_rank
    class(model_t), intent(in), target :: ml

    ! Broadcast the mode_t

    if (MPI_RANK /= root_rank) then
       md%ml => ml
    endif

    call bcast(md%mp, root_rank)
    call bcast(md%op, root_rank)

    if (MPI_RANK /= root_rank) then
       allocate(md%rt, SOURCE=c_rot_t(ml, md%mp, md%op))
    endif

    call bcast_alloc(md%x, root_rank)
    call bcast_alloc(md%y, root_rank)
    
    call bcast(md%x_ref, root_rank)
    call bcast(md%y_ref, root_rank)

    call bcast(md%omega, root_rank)

    call bcast(md%n_p, root_rank)
    call bcast(md%n_g, root_rank)
    call bcast(md%n_pg, root_rank)

    call bcast(md%n, root_rank)
    call bcast(md%n_iter, root_rank)

    call bcast(md%pruned, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

end module gyre_sol
