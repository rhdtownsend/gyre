! Module   : gyre_mode
! Purpose  : mode data
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_mode

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants
  use gyre_ext
  use gyre_freq
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory
  use gyre_sol
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure       :: ${NAME}_1_
    procedure       :: ${NAME}_v_
    generic, public :: ${NAME} => ${NAME}_1_, ${NAME}_v_
  $endsub

  type :: mode_t
     class(model_t), pointer     :: ml => null()
     type(sol_t), allocatable    :: sl
     class(c_rot_t), allocatable :: rt
     type(mode_par_t)            :: md_p
     type(osc_par_t)             :: os_p
     complex(WP)                 :: omega
     complex(WP)                 :: l_i
     complex(WP)                 :: scl
     integer, allocatable        :: s(:)
     real(WP), allocatable       :: x(:)
     integer                     :: s_ref
     real(WP)                    :: x_ref
     integer                     :: l
     integer                     :: m
     integer                     :: n_pg
     integer                     :: n_p
     integer                     :: n_g
     integer                     :: n_k
     logical                     :: pruned
   contains
     private
     procedure         :: op_assign_
     generic, public   :: assignment(=) => op_assign_
     procedure, public :: prune
     procedure, public :: freq
     $PROC_DECL(y)
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
     $PROC_DECL(lambda)
     $PROC_DECL(dE_dx)
     $PROC_DECL(dW_dx)
     $PROC_DECL(dC_dx)
     $PROC_DECL(F_j)
     $PROC_DECL(Yt_1)
     $PROC_DECL(Yt_2)
     $PROC_DECL(I_0)
     $PROC_DECL(I_1)
     $PROC_DECL(prop_type)
     procedure, public :: lag_T_eff
     procedure, public :: lag_g_eff
     procedure, public :: f_T
     procedure, public :: f_g
     procedure, public :: psi_T
     procedure, public :: psi_g
     procedure, public :: E
     procedure, public :: E_p
     procedure, public :: E_g
     procedure, public :: E_norm
     procedure, public :: W
     procedure, public :: C
     procedure, public :: omega_int
     procedure, public :: eta
     procedure         :: classify_
  end type mode_t

  ! Interfaces

  interface mode_t
     module procedure mode_t_
  end interface mode_t

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

  public :: mode_t
  public :: reallocate
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function mode_t_ (ml, sl, s, x, md_p, os_p) result (md)

    class(model_t), pointer, intent(in) :: ml
    type(sol_t), intent(in)             :: sl
    integer, intent(in)                 :: s(:)
    real(WP), intent(in)                :: x(:)
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(mode_t)                        :: md

    complex(WP) :: y_1_ref
    real(WP)    :: phase

    ! Construct the mode_t

    md%ml => ml

    allocate(md%sl, SOURCE=sl)
    allocate(md%rt, SOURCE=c_rot_t(ml, md_p, os_p))

    md%omega = sl%omega
    md%l_i = md%rt%l_i(md%omega)

    md%l = md_p%l
    md%m = md_p%m
    
    md%md_p = md_p
    md%os_p = os_p

    md%s = s
    md%x = x

    md%n_k = SIZE(s)

    md%pruned = .FALSE.

    ! Set the reference point

    call set_ref_()

    ! Normalize so that y_1 at the reference point is purely real, and
    ! the total mode energy is unity

    md%scl = 1._WP
    
    y_1_ref = md%y(1, md%s_ref, md%x_ref)

    phase = ATAN2(AIMAG(y_1_ref), REAL(y_1_ref))

    md%scl = 1._WP/SQRT(md%E())*EXP(CMPLX(0._WP, -phase, KIND=WP))

    ! Classify the mode

    call md%classify_()

    ! Finish

    return

  contains

    subroutine set_ref_ ()

      integer :: s

      ! Set up the reference point location

      md%x_ref = MIN(os_p%x_ref, ml%x_o(ml%n_s))

      seg_loop : do s = ml%n_s, 1, -1
         if (md%x_ref >= ml%x_i(s) .AND. md%x_ref <= ml%x_o(s)) then
            md%s_ref = s
            exit seg_loop
         end if
      end do seg_loop

      ! Finish

      return

    end subroutine set_ref_

  end function mode_t_

  !****

  subroutine op_assign_ (this, that)

    class(mode_t), intent(out) :: this
    type(mode_t), intent(in)   :: that

    ! Assign the mode_t (this routine shouldn't be necessary, but
    ! gfortran 4.9.x has memory issues using intrinsic assignment)

    this%ml => that%ml

    this%sl = that%sl
    allocate(this%rt, SOURCE=that%rt)

    this%omega = that%omega
    this%l_i = that%l_i

    this%l = that%l
    this%m = that%m
    
    this%md_p = that%md_p
    this%os_p = that%os_p

    this%s = that%s
    this%x = that%x

    this%n_k = that%n_k

    this%pruned = that%pruned

    this%s_ref = that%s_ref
    this%x_ref = that%x_ref

    this%scl = that%scl

    this%n_pg = that%n_pg
    this%n_p = that%n_p
    this%n_g = that%n_g

    ! Finish

    return

  end subroutine op_assign_

  !****

  $REALLOCATE(type(mode_t),1)

  !****

  $if($MPI)

  subroutine bcast_ (md, root_rank, ml)

    class(mode_t), intent(inout)       :: md
    integer, intent(in)                :: root_rank
    class(model_t), intent(in), target :: ml

    ! Broadcast the mode_t

    if (MPI_RANK /= root_rank) then
       md%ml => ml
    endif

    call bcast_alloc(md%sl, root_rank)

    call bcast(md%md_p, root_rank)
    call bcast(md%os_p, root_rank)

    if (MPI_RANK /= root_rank) then
       allocate(md%rt, SOURCE=c_rot_t(ml, md%md_p, md%os_p))
    endif

    call bcast(md%omega, root_rank)

    call bcast(md%l_i, root_rank)
    call bcast(md%scl, root_rank)

    call bcast_alloc(md%s, root_rank)
    call bcast_alloc(md%x, root_rank)
    call bcast(md%s_ref, root_rank)
    call bcast(md%x_ref, root_rank)

    call bcast(md%l, root_rank)
    call bcast(md%m, root_rank)

    call bcast(md%n_p, root_rank)
    call bcast(md%n_g, root_rank)
    call bcast(md%n_pg, root_rank)

    call bcast(md%n_k, root_rank)

    call bcast(md%pruned, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

  !****

  subroutine prune (this)

    class(mode_t), intent(inout) :: this

    ! Prune the mode_t

    if (.NOT. this%pruned) then

       deallocate(this%sl)

       deallocate(this%s)
       deallocate(this%x)

       this%pruned = .TRUE.

    endif

    ! Finish

    return

  end subroutine prune

  !****

  function freq (this, freq_units, freq_frame)

    class(mode_t), intent(in)          :: this
    character(*), intent(in)           :: freq_units
    character(*), optional, intent(in) :: freq_frame
    complex(WP)                        :: freq

    ! Calculate the frequency

    if (PRESENT(freq_frame)) then
       freq = freq_from_omega(this%omega, this%ml, freq_units, freq_frame, this%md_p, this%os_p)
    else
       freq = freq_from_omega(this%omega, this%ml, freq_units, 'INERTIAL', this%md_p, this%os_p)
    endif

    ! Finish
    
    return

  end function freq

  !****

  function y_1_ (this, i, s, x) result (y)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: i
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: y

    $ASSERT(.NOT. this%pruned,Cannot evaluate y from pruned data)

    ! Evaluate y(i)

    y = this%scl*this%sl%y(i, s, x)

    ! Finish

    return

  end function y_1_

  !****

  function xi_r_1_ (this, s, x) result (xi_r)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: xi_r

    complex(WP) :: y_1

    ! Evaluate the radial displacement perturbation, in units of
    ! R_star

    associate (l_i => this%l_i)

      y_1 = this%y(1, s, x)

      if (l_i /= 1._WP) then

         if (x /= 0._WP) then
            xi_r = y_1*x**(l_i-1._WP)
         else
            xi_r = 0._WP
         endif
       
      else
       
         xi_r = y_1

      endif

    end associate

    ! Finish

    return

  end function xi_r_1_

  !****

  function xi_h_1_ (this, s, x) result (xi_h)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: xi_h

    complex(WP) :: y_2
    complex(WP) :: y_3
    real(WP)    :: c_1
    complex(WP) :: omega_c

    ! Evaluate the horizontal displacement perturbation, in units of
    ! R_star

    if (this%l /= 0) then

       associate (l_i => this%l_i)

         y_2 = this%y(2, s, x)
         y_3 = this%y(3, s, x)

         c_1 = this%ml%c_1(s, x)

         omega_c = this%rt%omega_c(s, x, this%omega)
      
         if (l_i /= 1._WP) then

            if (x /= 0._WP) then
               xi_h = (y_2+y_3)*x**(l_i-1._WP)/(c_1*omega_c**2)
            else
               xi_h = 0._WP
            end if

         else
            
            xi_h = (y_2+y_3)/(c_1*omega_c**2)

         endif

       end associate

    else

       xi_h = 0._WP

    end if

    ! Finish

    return
    
  end function xi_h_1_

  !****

  function eul_phi_1_ (this, s, x) result (eul_phi)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: eul_phi

    complex(WP) :: y_3
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian gravitational potential perturbation, in
    ! units of G M_star / R_star

    associate (l_i => this%l_i)

      y_3 = this%y(3, s, x)

      c_1 = this%ml%c_1(s, x)

      if (l_i /= 0._WP) then

         if (x /= 0._WP) then
            eul_phi = y_3*x**l_i/c_1
         else
            eul_phi = 0._WP
         endif

      else

         eul_phi = y_3/c_1

      endif

    end associate

    ! Finish

    return

  end function eul_phi_1_

  !****

  function deul_phi_1_ (this, s, x) result (deul_phi)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: deul_phi

    complex(WP) :: y_4
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian potential gradient (gravity) perturbation,
    ! in units of G M_star / R_star**2

    associate (l_i => this%l_i)

      y_4 = this%y(4, s, x)

      c_1 = this%ml%c_1(s, x)

      if (l_i /= 1._WP) then

         if (x /= 0._WP) then
            deul_phi = y_4*x**(l_i-1._WP)/c_1
         else
            deul_phi = 0._WP
         end if
       
      else

         deul_phi = y_4/c_1

      end if

    end associate

    ! Finish

    return

  end function deul_phi_1_

  !****

  function lag_S_1_ (this, s, x) result (lag_S)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lag_S

    complex(WP) :: y_5

    ! Evaluate the Lagrangian specific entropy perturbation, in units
    ! of c_p

    associate (l_i => this%l_i)

      y_5 = this%y(5, s, x)

      if (x /= 0._WP) then
         lag_S = y_5*x**(l_i-2._WP)
      else
         lag_S = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_S_1_

  !****

  function lag_L_1_ (this, s, x) result (lag_L)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lag_L

    complex(WP) :: y_6

    ! Evaluate the Lagrangian radiative luminosity perturbation, in
    ! units of L_star

    associate (l_i => this%l_i)

      y_6 = this%y(6, s, x)

      if (x /= 0._WP) then
         lag_L = y_6*x**(l_i+1._WP)
      else
         lag_L = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_L_1_

  !****

  function eul_P_1_ (this, s, x) result (eul_P)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: eul_P

    complex(WP) :: xi_r
    complex(WP) :: lag_P
    real(WP)    :: V_2

    ! Evaluate the Eulerian pressure perturbation, in units of P

    xi_r = this%xi_r(s, x)
    lag_P = this%lag_P(s, x)

    V_2 = this%ml%V_2(s, x)

    eul_P = lag_P + V_2*x*xi_r

    ! Finish

    return

  end function eul_P_1_
  
  !****

  function lag_P_1_ (this, s, x) result (lag_P)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lag_P

    complex(WP) :: y_1
    complex(WP) :: y_2
    real(WP)    :: V_2

    ! Evaluate the Lagrangian pressure perturbation, in units of P

    associate (l_i => this%l_i)

      y_1 = this%y(1, s, x)
      y_2 = this%y(2, s, x)

      V_2 = this%ml%V_2(s, x)

      if (l_i /= 0._WP) then

         if (x /= 0._WP) then
            lag_P = V_2*(y_2 - y_1)*x**l_i
         else
            lag_P = 0._WP
         end if

      else

         lag_P = V_2*(y_2 - y_1)

      endif

    end associate

    ! Finish

    return

  end function lag_P_1_

  !****

  function eul_rho_1_ (this, s, x) result (eul_rho)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: eul_rho

    complex(WP) :: xi_r
    complex(WP) :: lag_rho
    real(WP)    :: U
    real(WP)    :: dU
    real(WP)    :: D

    ! Evaluate the Eulerian density perturbation, in units of rho

    xi_r = this%xi_r(s, x)
    lag_rho = this%lag_rho(s, x)

    U = this%ml%U(s, x)
    dU = this%ml%dU(s, x)

    D = dU + U - 3._WP

    if (x /= 0._WP) then
       eul_rho = lag_rho - D*xi_r/x
    else
       eul_rho = lag_rho
    endif

    ! Finish

    return

  end function eul_rho_1_

  !****

  function lag_rho_1_ (this, s, x) result (lag_rho)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lag_rho

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: Gamma_1
    real(WP)    :: delta

    ! Evaluate the Lagrangian density perturbation, in units of
    ! rho. This expression implements eqn. 13.83 of [Unn1989]

    lag_P = this%lag_P(s, x)
    lag_S = this%lag_S(s, x)

    Gamma_1 = this%ml%Gamma_1(s, x)
    delta = this%ml%delta(s, x)

    lag_rho = lag_P/Gamma_1 - delta*lag_S

    ! Finish

    return

  end function lag_rho_1_

  !****

  function eul_T_1_ (this, s, x) result (eul_T)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: eul_T

    complex(WP) :: xi_r
    complex(WP) :: lag_T
    real(WP)    :: V_2
    real(WP)    :: nabla

    ! Evaluate the Lagrangian temperature perturbation, in units of T

    xi_r = this%xi_r(s, x)
    lag_T = this%lag_T(s, x)

    V_2 = this%ml%V_2(s, x)
    nabla = this%ml%nabla(s, x)
      
    eul_T = lag_T + nabla*V_2*x*xi_r

    ! Finish

    return

  end function eul_T_1_

  !****

  function lag_T_1_ (this, s, x) result (lag_T)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lag_T

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: nabla_ad

    ! Evaluate the Lagrangian temperature perturbation, in units of
    ! T. This expression implements eqn. 13.84 of [Unn1989]

    lag_P = this%lag_P(s, x)
    lag_S = this%lag_S(s, x)

    nabla_ad = this%ml%nabla_ad(s, x)
      
    lag_T = nabla_ad*lag_P + lag_S

    ! Finish

    return

  end function lag_T_1_

  !****

  function lambda_1_ (this, s, x) result (lambda)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: lambda

    ! Evaluate the angular eigenvalue

    lambda = this%rt%lambda(s, x, this%omega)

    ! Finish

    return

  end function lambda_1_
    
  !****

  function dE_dx_1_ (this, s, x) result (dE_dx)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    real(WP)                  :: dE_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: U
    real(WP)    :: c_1

    ! Evaluate the differential mode inertia, in units of M_star
    ! R_star**2. This expression is based on eqn. 3.139 of [Aer2010],
    ! with the initial factor of 4 pi cancelled by their definitions of
    ! \tilde{\xi}_r and \tilde{\xi}_h (cf. eqns. 3.124 and 3.131, ibid)

    xi_r = this%xi_r(s, x)
    xi_h = this%xi_h(s, x)

    lambda = this%lambda(s, x)

    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    dE_dx = (ABS(xi_r)**2 + ABS(lambda)*ABS(xi_h)**2)*U*x**2/(4._WP*PI*c_1)

    ! Finish

    return

  end function dE_dx_1_

  !****

  function dW_dx_1_ (this, s, x) result (dW_dx)

    use gyre_evol_model

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    real(WP)                  :: dW_dx

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

    lag_T = this%lag_T(s, x)
    lag_S = this%lag_S(s, x)
    
    c_thm = this%ml%c_thm(s, x)

    dW_dx = PI*AIMAG(CONJG(lag_T)*lag_S)*c_thm*x**2*t_dyn/t_kh

    ! Finish

    return

  end function dW_dx_1_

  !****

  function dW_dx_eps_1_ (this, s, x) result (dW_dx_eps)

    use gyre_evol_model

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    real(WP)                  :: dW_dx_eps

    ! Evaluate the differential work associated with nuclear
    ! processes, in units of G M_star**2/R_star.  This expression is
    ! based on eqn. X.XX of [Unn1989]

    dW_dx_eps = 0._WP

    $ABORT(Not yet implemented)

    ! Finish

    return

  end function dW_dx_eps_1_

  !****

  function dC_dx_1_ (this, s, x) result (dC_dx)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    real(WP)                  :: dC_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: lambda
    real(WP)    :: E
    
    ! Calculate the (unnormalized) rotation splitting kernel. This is
    ! based on based on equations 3.356 & 3.357 of [Aer2010], with
    ! dC_dx = beta*K

    associate (l => this%l)

      xi_r = this%xi_r(s, x)
      xi_h = this%xi_h(s, x)

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      lambda = this%lambda(s, x)

      ! (The following holds due to the standard normalization)
 
      E = 1._WP

      dC_dx = REAL((ABS(xi_r)**2 + (lambda-1._WP)*ABS(xi_h)**2 - &
                   2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1)/E

    end associate

    ! Finish

    return

  end function dC_dx_1_

  !****

  function F_j_1_ (this, s, x) result (F_j)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    real(WP)                  :: F_j

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: c_1
    real(WP)    :: U
    complex(WP) :: omega_c
    
    ! Evaluate the angle-averaged angular momentum flux due to
    ! Reynolds stress, in units of G M_star**2/R_star**3.  This
    ! expression is based on eqn. 21 of [LeeSai1993]

    associate (m => this%m)

      xi_r = this%xi_r(s, x)
      xi_h = this%xi_h(s, x)

      c_1 = this%ml%c_1(s, x)
      U = this%ml%U(s, x)

      omega_c = this%rt%omega_c(s, x, this%omega)

      F_j = -m*ABS(omega_c**2)*x*U*AIMAG(CONJG(xi_r)*xi_h)/(32._WP*PI**2*c_1)

    end associate

    ! Finish

    return

  end function F_j_1_
  
  !****

  function Yt_1_1_ (this, s, x) result (Yt_1)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: Yt_1

    complex(WP) :: y_1
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: J

    ! Evaluate the Takata Y_1 function. This expression is equivalent to
    ! eqn. 69 of [Tak2006b], divided by x**(2-l)

    y_1 = this%y(1, s, x)
    y_3 = this%y(3, s, x)
    y_4 = this%y(4, s, x)

    J = 1._WP - this%ml%U(s, x)/3._WP

    Yt_1 = J*y_1 + (y_3 - y_4)/3._WP

    ! Finish

    return

  end function Yt_1_1_

  !****

  function Yt_2_1_ (this, s, x) result (Yt_2)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: Yt_2

    complex(WP) :: y_1
    complex(WP) :: y_2

    ! Evaluate the Takata Y_2 function. This expression is equivalent to 
    ! eqn. 70 of [Tak2006b], divided by V

    y_1 = this%y(1, s, x)
    y_2 = this%y(2, s, x)

    Yt_2 = y_2 - y_1

    ! Finish

    return

  end function Yt_2_1_

  !****

  function I_0_1_ (this, s, x) result (I_0)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: I_0

    complex(WP) :: y_1
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1
    
    ! Evaluate the I_0 integral, which should be zero for radial
    ! modes. This expression is based on eqn. 42 of [Tak2006a]

    associate (l_i => this%l_i)

      y_1 = this%y(1, s, x)
      y_4 = this%y(4, s, x)

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      if (x /= 0._WP) then
         I_0 = x**(l_i+1._WP)*(U*y_1 + y_4)/c_1
      else
         I_0 = 0._WP
      endif

    end associate

    ! Finish

    return

  end function I_0_1_

  !****

  function I_1_1_ (this, s, x) result (I_1)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    complex(WP)               :: I_1

    complex(WP) :: y_1
    complex(WP) :: y_2
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: omega_c

    ! Evaluate the I_0 integral, which should be zero for dipole
    ! modes. This expression is based on eqn. 43 of [Tak2006a]

    associate (l_i => this%l_i)

      y_1 = this%y(1, s, x)
      y_2 = this%y(2, s, x)
      y_3 = this%y(3, s, x)
      y_4 = this%y(4, s, x)

      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      omega_c = this%rt%omega_c(s, x, this%omega)

      if (x /= 0._WP) then
         I_1 = x**(l_i+2._WP)*(c_1*omega_c**2*U*y_1 - U*y_2 + &
               (c_1*omega_c**2 - 2._WP)*y_3 + (c_1*omega_c**2 - 1._WP)*y_4)/c_1**2
      else
         I_1 = 0._WP
      endif

    end associate
      
    ! Finish

    return

  end function I_1_1_

  !****

  function prop_type_1_ (this, s, x) result (prop_type)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s
    real(WP), intent(in)      :: x
    integer                   :: prop_type

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0
    real(WP) :: gamma

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    if (this%ml%vacuum(s, x)) then

       prop_type = 0

    else

       ! Calculate the discriminant gamma

       V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
       As = this%ml%As(s, x)
       U = this%ml%U(s, x)
       c_1 = this%ml%c_1(s, x)

       lambda = REAL(this%lambda(s, x))

       omega_c = REAL(this%rt%omega_c(s, x, this%omega))

       g_4 = -4._WP*V_g*c_1
       g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
       g_0 = -4._WP*lambda*As/c_1
       
       gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

       ! Use the sign of gamma to set up prop_type

       if (gamma > 0._WP) then

          prop_type = 0

       else

          if (REAL(c_1*omega_c**2) < As) then
             prop_type = -1
          else
             prop_type = 1
          endif

       endif

    endif
    
    ! Finish

    return

  end function prop_type_1_

  !****

  function y_v_ (this, i, s, x) result (y)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: i
    integer, intent(in)       :: s(:)
    real(WP), intent(in)      :: x(:)
    complex(WP)               :: y(SIZE(s))

     integer :: k

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Evaluate y(i)

    !$OMP PARALLEL DO
    do k = 1, SIZE(s)
       y(k) = this%y(i, s(k), x(k))
    end do

    ! Finish

    return

  end function y_v_

  !****

  $define $PROC_V $sub

  $local $NAME $1
  $local $TYPE $2

  function ${NAME}_v_ (this, s, x) result (${NAME})

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: s(:)
    real(WP), intent(in)      :: x(:)
    $TYPE                     :: $NAME(SIZE(s))

     integer :: k

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Evaluate $NAME

    !$OMP PARALLEL DO
    do k = 1, SIZE(s)
       ${NAME}(k) = this%${NAME}(s(k), x(k))
    end do

    ! Finish

    return

  end function ${NAME}_v_

  $endsub

  $PROC_V(xi_r,complex(WP))
  $PROC_V(xi_h,complex(WP))
  $PROC_V(eul_phi,complex(WP))
  $PROC_V(deul_phi,complex(WP))
  $PROC_V(lag_S,complex(WP))
  $PROC_V(lag_L,complex(WP))
  $PROC_V(eul_P,complex(WP))
  $PROC_V(lag_P,complex(WP))
  $PROC_V(eul_rho,complex(WP))
  $PROC_V(lag_rho,complex(WP))
  $PROC_V(eul_T,complex(WP))
  $PROC_V(lag_T,complex(WP))
  $PROC_V(lambda,complex(WP))
  $PROC_V(dE_dx,real(WP))
  $PROC_V(dW_dx,real(WP))
  $PROC_V(dC_dx,real(WP))
  $PROC_V(F_j,real(WP))
  $PROC_V(Yt_1,complex(WP))
  $PROC_V(Yt_2,complex(WP))
  $PROC_V(I_0,complex(WP))
  $PROC_V(I_1,complex(WP))
  $PROC_V(prop_type,integer)

  !****

  function lag_T_eff (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_T_eff

    complex(WP) :: xi_r
    complex(WP) :: lag_L

    ! Calculate the effective temperature perturbation at x_ref
    ! (assumed to correspond to the photosphere), in units of
    ! T_eff. This expression is based on the standard definition of
    ! effective temperature

    associate (s => this%s_ref, x => this%x_ref)

      xi_r = this%xi_r(s, x)
      lag_L = this%lag_L(s, x)

      lag_T_eff = 0.25_WP*(lag_L - 2._WP*xi_r)

    end associate

    ! Finish

    return

  end function lag_T_eff

  !****

  function lag_g_eff (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_g_eff

    complex(WP) :: xi_r
    complex(WP) :: deul_phi
    real(WP)    :: c_1
    real(WP)    :: U

    ! Calculate the effective gravity perturbation at x_ref (assumed
    ! to correspond to the photosphere), in units of the gravity. This
    ! expression is based on eqn. 24 of [Dup2002]

    associate (s => this%s_ref, x => this%x_ref, omega => this%omega)

      xi_r = this%xi_r(s, x)
      deul_phi = this%deul_phi(s, x)

      c_1 = this%ml%c_1(s, x)
      U = this%ml%c_1(s, x)

      lag_g_eff = (c_1/x)*deul_phi + (U - (2._WP + c_1*omega**2))*xi_r/x

    end associate

    ! Finish

    return

  end function lag_g_eff

  !****

  function f_T (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_T

    complex(WP) :: C_T

    ! Evaluate the non-adiabatic f_T parameter. This is expression is
    ! based on eqn. 5 of [Dup2003]

    associate (s => this%s_ref, x => this%x_ref)

      C_T = this%lag_T_eff()/this%xi_r(s, x)

      f_T = ABS(C_T)

    end associate

    ! Finish

    return

  end function f_T

  !****

  function f_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_g

    complex(WP) :: C_g

    ! Evaluate the non-adiabatic f_g parameter. This is expression is
    ! based on eqn. 6 of [Dup2003]

    associate (s => this%s_ref, x => this%x_ref)

      C_g = this%lag_g_eff()/this%xi_r(s, x)

      f_g = -ABS(C_g)

    end associate

    ! Finish

    return

  end function f_g

  !****

  function psi_T (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: psi_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic psi_T parameter, in radians. This is
    ! expression is based on eqn. 5 of [Dup2003]

    associate (s => this%s_ref, x => this%x_ref)

      C_T = this%lag_T_eff()/this%xi_r(s, x)

      psi_T = ATAN2(AIMAG(C_T), REAL(C_T))

    end associate

    ! Finish

    return

  end function psi_T

  !****

  function psi_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: psi_g

    ! Calculate the non-adiabatic psi_g parameter, in radians

    psi_g = PI

    ! Finish

    return

  end function psi_g

  !****

  function E (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E
    
    ! Calculate the mode inertia, in units of M_star R_star**2

    associate (s => this%s, x => this%x)

      E = integrate(x, this%dE_dx(s, x))

    end associate

    ! Finish

    return

  end function E

  !****

  function E_p (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E_p

    ! Calculate the mode inertia in acoustic-wave propagation regions,
    ! in units of M_star R_star**2

    associate (s => this%s, x => this%x)

      E_p = integrate(x, this%dE_dx(s, x), mask=(this%prop_type(s, x) == 1))

    end associate
    
    ! Finish

    return

  end function E_p

  !****

  function E_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E_g

    ! Calculate the mode inertia in gravity-wave propagation regions,
    ! in units of M_star R_star**2

    associate (s => this%s, x => this%x)

      E_g = integrate(x, this%dE_dx(s, x), mask=(this%prop_type(s, x) == -1))

    end associate
    
    ! Finish

    return

  end function E_g

  !****

  function E_norm (this)
 
    class(mode_t), intent(in) :: this
    real(WP)                  :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: A2

    ! Calculate the normalized mode inertia. This expression is based
    ! on eqn. 3.140 of [Aer2010]

    associate (s => this%s_ref, x => this%x_ref)

      ! (The following holds due to the standard normalization)

      E = 1._WP

      xi_r = this%xi_r(s, x)
      xi_h = this%xi_h(s, x)

      lambda = this%lambda(s, x)

      select case (this%os_p%inertia_norm)
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

    end associate

    ! Finish

    return

  end function E_norm

  !****

  function W (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: W
    
    ! Calculate the total work, in units of G M_star**2/R_star

    associate (s => this%s, x => this%x)

      W = integrate(x, this%dW_dx(s, x))

    end associate

    ! Finish

    return

  end function W

  !****

  function C (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: C
    
    ! Calculate the Ledoux rotational splitting coefficient

    associate (s => this%s, x => this%x)

      C = integrate(x, this%dC_dx(s, x))

    end associate

    ! Finish

    return

  end function C

  !****

  function omega_int (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: omega_int

    integer     :: k
    complex(WP) :: xi_r
    complex(WP) :: eul_phi
    complex(WP) :: eul_rho
    complex(WP) :: lag_rho
    complex(WP) :: lag_P
    real(WP)    :: V_2
    real(WP)    :: As
    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: Gamma_1
    real(WP)    :: V_g
    real(WP)    :: x4_V
    complex(WP) :: f_th(this%n_k)
    complex(WP) :: f_re(this%n_k)
    complex(WP) :: f_gr(this%n_k)
    complex(WP) :: f_xi(this%n_k)
    complex(WP) :: W_th
    complex(WP) :: W_re
    complex(WP) :: W_gr
    complex(WP) :: W_xi

    ! Calculate the dimensionless frequency from the integral
    ! expression in eqn. (1.71) of [Dup2003]

    do k = 1, this%n_k

       associate (s => this%s(k), x => this%x(k))

         xi_r = this%xi_r(s, x)
         eul_phi = this%eul_phi(s, x)
         eul_rho = this%eul_rho(s, x)
         lag_rho = this%lag_rho(s, x)
         lag_P = this%lag_P(s, x)

         V_2 = this%ml%V_2(s, x)
         As = this%ml%As(s, x)
         U = this%ml%U(s, x)
         c_1 = this%ml%c_1(s, x)

         Gamma_1 = this%ml%Gamma_1(s, x)

         V_g = V_2*x**2/Gamma_1
         x4_V = x**2/V_2

         f_th(k) = CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2))
         f_re(k) = 2._WP*REAL(lag_rho*CONJG(xi_r)*(x/c_1)*(x**2*U/c_1))
         f_gr(k) = CONJG(eul_rho)*eul_phi*(x**2*U/c_1)
         f_xi(k) = -ABS(xi_r)**2*(x/c_1)*(x*U*(-V_g-As)/c_1)

       end associate

    end do

    W_th = integrate(this%x, f_th)
    W_re = integrate(this%x, f_re)
    W_gr = integrate(this%x, f_gr)
    W_xi = integrate(this%x, f_xi)

    omega_int = SQRT(4._WP*PI*(W_th + W_re + W_gr + W_xi)/this%E())

    ! Finish

    return

  end function omega_int

  !****

  function eta (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: eta

    real(WP) :: dW_dx(this%n_k)

    ! Calculate the normalized growth rate defined by [Stel1978]

    associate (s => this%s, x => this%x)

      dW_dx = this%dW_dx(s, x)

      eta = integrate(x, dW_dx)/integrate(x, ABS(dW_dx))

    end associate

    ! Finish

    return

  end function eta

  !****

  subroutine classify_ (this)

    class(mode_t), intent(inout) :: this

    real(WP) :: y_1(this%n_k)
    real(WP) :: y_2(this%n_k)
    integer  :: k_i
    integer  :: k_o
    integer  :: n_c
    integer  :: n_a

    ! Classify the mode based on its eigenfunctions

    if (this%l == 0) then

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       y_1 = REAL(this%y(1, this%s, this%x))
       y_2 = REAL(this%y(2, this%s, this%x))

       mono_loop : do k_i = 2, this%n_k-1
          if ((y_1(k_i) >= y_1(k_i-1) .AND. y_1(k_i+1) >= y_1(k_i)) .OR. &
              (y_1(k_i) <= y_1(k_i-1) .AND. y_1(k_i+1) <= y_1(k_i))) exit mono_loop
       end do mono_loop

       ! Count winding numbers

       call count_windings_(y_1(k_i:), y_2(k_i:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       this%n_p = n_a + 1
       this%n_g = n_c

       this%n_pg = this%n_p - this%n_g

    elseif (this%l == 1 .AND. .NOT. this%os_p%cowling_approx) then

       ! Dipole modes (non-Cowling)

       ! Set up the Takata Y^a_1 and Y^a_2 functions

       y_1 = REAL(this%Yt_1(this%s, this%x))
       y_2 = REAL(this%Yt_2(this%s, this%x))

       ! Count winding numbers, taking care to avoid counting nodes at
       ! the center and surface (check x(1) rather than Yt_1(1),
       ! because the latter can be non-zero at x=0 due to rounding
       ! errors)

       if (this%x(1) == 0._WP) then
          k_i = 2
       else
          k_i = 1
       endif

       if (y_1(this%n_k) == 0._WP) then
          k_o = this%n_k-1
       else
          k_o = this%n_k
       endif
       
       call count_windings_(y_1(k_i:k_o), y_2(k_i:k_o), n_c, n_a)

       ! Classify

       this%n_p = n_a
       this%n_g = n_c

       if (this%n_p >= this%n_g) then
          this%n_pg = this%n_p - this%n_g + 1
       else
          this%n_pg = this%n_p - this%n_g
       endif

    else

       ! Other modes

       y_1 = REAL(this%y(1, this%s, this%x))
       y_2 = REAL(this%y(2, this%s, this%x) + this%y(3, this%s, this%x))

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

end module gyre_mode
