  ! Module   : core_multi_func (version 1.0)
! Purpose  : systems of functions, evaluation and root
!          : currently only handles real numbers (no complex)
!          : solves R^n => R^n systems
!          : Based on J.E. Dennis, R.B. Schnabel, "Numerical Methods for.."
!             (1983/1996)
!          : and Numerical Recipes (Press et al)

$include 'core.inc'

module core_multi_func

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Parameters

  real(WP), parameter :: TOL = 1e-9_WP
  integer, parameter  :: MAX_ITER = 800
  real(WP), parameter :: MAXSTEP_FAC = 1e3_WP
  real(WP), parameter :: STEP_TOL = 1e-12_WP
  real(WP), parameter :: NORM_TOL = 1e-12_WP

  ! termination codes meaning
  integer, parameter :: ALL_OK = 1
  integer, parameter :: TOO_SLOW = 2
  integer, parameter :: NORM_DIDNT_DECREASE = 3
  integer, parameter :: MAY_NEED_DIFFERENT_INITIAL = 6
  integer, parameter :: MAX_ITER_EXCEEDED = 4
  integer, parameter :: QUITE_LARGE_STEPS = 5
  integer, parameter :: EVAL_RETURNED_BAD = 7

  ! global convergence options
  integer, parameter :: LINE_SEARCH = 1
  integer, parameter :: HOOK_SEARCH = 2
  integer, parameter :: DOGLEG_SEARCH = 3

  ! evaluation status
  integer, parameter :: EVAL_OK  = 0
  integer, parameter :: EVAL_BAD = 1
  
  ! Derived-type definitions

  type, abstract :: multi_func_t
     integer                      :: neq ! number of equations and unknowns
     logical                      :: broyden ! use secant method to estimate jac
     integer                      :: global ! options for global convergence
   contains
     private
     procedure(eval_r_), deferred :: eval_r_
     generic, public              :: eval => eval_r_
     procedure, public            :: root => root_r_
     procedure                    :: jac=>jac_r_
     procedure                    :: gradient => gradient_
     procedure                    :: eval_norm => eval_norm_
     procedure                    :: line_search => line_search_
     procedure                    :: hook_search => hook_search_
     procedure                    :: hook_step => hook_step_
     procedure                    :: nonlinear_solve => nonlinear_solve_
     procedure                    :: nonlinear_solve_unfac => nonlinear_solve_unfac_
     procedure                    :: trust_region_update => trust_region_update_
     procedure                    :: QR_solve => QR_solve_
     procedure                    :: QR_decomp => QR_decomp_
     procedure                    :: QR_update => QR_update_
     procedure                    :: Q_form => Q_form_
     procedure                    :: condest => condest_
     procedure                    :: rsolve => rsolve_
     procedure                    :: cholsolve => cholsolve_
     procedure                    :: choldecomp_pd => choldecomp_pd_
     procedure                    :: lsolve => lsolve_
     procedure                    :: ltsolve => ltsolve_
     procedure                    :: broyden_fac => broyden_fac_
     procedure                    :: jac_rotate => jac_rotate_
     procedure                    :: stop_qq => stop_qq_
  end type multi_func_t

  ! Interfaces

  abstract interface
    subroutine eval_r_ (this, x, F_x, status)
      use core_kinds
      import multi_func_t
      class(multi_func_t), intent(inout) :: this
      real(WP), intent(in)           :: x(:)
      real(WP), intent(inout)        :: F_x(:)
      integer, intent(out), optional :: status 
    end subroutine eval_r_
  end interface

  ! Access specifiers

  private

  public :: multi_func_t
  public :: LINE_SEARCH
  public :: HOOK_SEARCH
  public :: DOGLEG_SEARCH
  public :: EVAL_BAD
  public :: EVAL_OK
  ! Procedures

contains

  ! based on algorithm A6.5.1fac
  
  ! find root of system given initial guess
  ! Make sure initial guess is good!
  
  function root_r_ (this, x_init, term_code, x_typical, F_typical, n_iter, broyden, global) result (x_out)

    class(multi_func_t), intent(inout)    :: this
    real(WP), intent(in)                  :: x_init(:)
    integer, intent(out)                  :: term_code
    real(WP), optional, intent(in)        :: x_typical(:)
    real(WP), optional, intent(in)        :: F_typical(:)
    integer, optional, intent(in)         :: n_iter
    logical, optional, intent(in)         :: broyden
    integer, optional, intent(in)         :: global
    
    real(WP), allocatable                 :: x(:)
    real(WP), allocatable                 :: x_out(:)
    real(WP), allocatable                 :: x_n(:)
    real(WP), allocatable                 :: F_a(:)
    real(WP), allocatable                 :: F_n(:)
    real(WP)                              :: f_norm_a
    real(WP)                              :: f_norm_n

    real(WP), allocatable                 :: g_a(:)
    real(WP), allocatable                 :: g_temp(:)
    real(WP), allocatable                 :: M(:,:)
    real(WP), allocatable                 :: M2(:)
    real(WP), allocatable                 :: J_a(:,:)
    real(WP), allocatable                 :: H(:,:)
    real(WP), allocatable                 :: s_n(:)

    real(WP), allocatable                 :: x_scale_(:)
    real(WP), allocatable                 :: F_scale_(:)
    
    real(WP) :: check
    real(WP) :: fafs
    real(WP) :: delta
    real(WP) :: delta_prev
    real(WP) :: mu
    real(WP) :: phi
    real(WP) :: phi_prime
    
    integer :: nx
    integer :: nsc
    logical :: max_taken
    integer :: consec_max
    integer :: i
    integer :: j
    integer :: k
    integer :: iter_count
    integer :: itermax

    integer :: ret_status
    integer :: eval_status
    
    logical :: no_sol
    logical :: restart
    logical :: to_secant
    
    nx = SIZE(x_init)
    $ASSERT(nx == this%neq, Wrong size of initial guess array)

    ! allocate needed arrays
       
    allocate(x_scale_(this%neq))
    allocate(F_scale_(this%neq)) 
    allocate(x(nx))
    allocate(x_out(nx))
    allocate(x_n(nx))
    allocate(F_a(nx))
    allocate(F_n(nx))
    allocate(J_a(nx,nx))
    allocate(H(nx,nx))
    allocate(s_n(nx))
    allocate(g_a(nx))
    allocate(g_temp(nx))
    allocate(M(nx,nx))
    allocate(M2(nx))

    if(PRESENT(broyden)) then
       this%broyden = broyden
    else
       this%broyden = .TRUE. ! use secant method by default
    end if
    
    if(PRESENT(global)) then
       $ASSERT(global > 0 .AND. global < 4, wrong global convergence option)
       this%global = global
    else
       this%global = LINE_SEARCH ! use line_search method by default
    end if
    
    ! _a is current position, _n is new position
    
    call this%eval_r_(x_init, F_a, eval_status)
    x = x_init

    if(eval_status == EVAL_BAD) then
       print *, 'init guess bad'
       term_code = EVAL_RETURNED_BAD
       x_out = x
       return
    end if
       
    
    ! check for root

    check = 0._WP
    do i=1,nx
       check = MAX(ABS(F_a(i)), check)
    end do

    if(check < TOL) then
       ! root is at x_init
       return
    end if

    delta = -1._WP
    
    ! else actually solve for root

    ! first check optional scaling matrices

    ! typical(i) is the typical size of x(i)
    ! scaling helps algorithm if x(i) is several orders
    !  of magnitude different than x(j) (better to be similar mag)
    
    if(PRESENT(x_typical)) then
       nsc = SIZE(x_typical)
       $ASSERT(nsc == this%neq, Wrong size of x typical array)

       x_scale_ = 1._WP/ABS(x_typical)
    else
       x_scale_ = 1._WP
    end if

    if(PRESENT(F_typical)) then
       nsc = SIZE(f_typical)
       $ASSERT(nsc == this%neq, Wrong size of f typical array)

       F_scale_ = 1._WP/ABS(f_typical)
    else
       F_scale_ = 1._WP
    end if

    ! set up initial 

    
    call this%jac(x, F_a, x_scale_, J_a)

    f_norm_a = 0.5*this%eval_norm(F_a, F_scale_)**2
    
    g_a = 0._WP
    
    do j=1,nx
       fafs = F_a(j)*F_scale_(j)**2
       do i=1,nx
          g_a(i) = g_a(i) + J_a(j,i)*fafs
       end do
    end do

    if(PRESENT(n_iter)) then
       itermax = n_iter
    else
       itermax = MAX_ITER
    end if

    ! iterate

    term_code = 0
    iter_count = 0
    ret_status = 0
    restart = .TRUE.
    to_secant = .FALSE.
    no_sol = .TRUE.
    
    do while(no_sol)
       if(.NOT.(to_secant)) then
          !print *, 'i', iter_count
          iter_count = iter_count + 1

          if(this%broyden) then
             call this%nonlinear_solve(restart, F_a, g_a, M, M2, J_a, H, s_n, F_scale_, x_scale_)
          else
             call this%nonlinear_solve_unfac(F_a, g_a, J_a, M, H, s_n, F_scale_, x_scale_)
          end if

          ! outputs x_n, F_n, f_norm_n

          if(this%global == LINE_SEARCH) then
             call this%line_search(x, x_n, f_norm_a, f_norm_n, F_n, g_a, s_n, ret_status, max_taken, F_scale_, x_scale_)
          elseif(this%global == HOOK_SEARCH) then
             call this%hook_search(x, x_n, f_norm_a, f_norm_n, F_n, g_a, s_n, M, H, F_scale_, x_scale_, iter_count, delta, delta_prev, mu, phi, phi_prime, ret_status, max_taken)
          elseif(this%global == DOGLEG_SEARCH) then
             $ABORT(dogleg search not implemented)
          else
             $ABORT(something totally crapped out here)
          end if

          !print *, 'ls status', ret_status, restart
          !print *,'old', x, F_a
          !print *,'new', x_n, F_n
          
          if(max_taken) then
             consec_max = consec_max + 1
          else
             consec_max = 0
          end if
          
       end if

       if(this%broyden .AND. ((ret_status == 1 .AND. .NOT.(restart)) .OR. to_secant)) then

          !print *, 'new jac'
          
          ! status-> need to restart secant method
          call this%jac(x, F_a, x_scale_, J_a)

          g_a = 0._WP
    
          do j=1,nx
             fafs = F_a(j)*F_scale_(j)**2
             do i=1,nx
                g_a(i) = g_a(i) + J_a(j,i)*fafs
             end do
          end do

          to_secant = .FALSE.
          restart = .TRUE.

          if(this%global > LINE_SEARCH) then
             delta = -1._WP
          end if

          !skip rest of do while

       else
          ! status-> continue
          ! update jacobian and related g-array

          if(this%broyden) then
          
             if(SUM(x_n-x) == 0._WP) then
                ! do nothing; program will terminate with code 3
                !  for failure to find smaller norm
             else
                call this%broyden_fac(x, x_n, F_a, F_n, J_a, M, M2, F_scale_, x_scale_)
                ! calculate g using QR factoriaztion
                
                g_a = 0._WP
                do j=1,nx
                   fafs = F_n(j)*F_scale_(j)
                   do i=1,nx
                      g_a(i) = g_a(i) + fafs*J_a(i,j)
                   end do
                end do
                
                ! QR part 2
                
                g_temp = g_a
                g_a = 0._WP
                
                do i=nx,1,-1
                   do j=1,i
                      g_a(i) = g_a(i) + M(j,i)*g_temp(j) !sdasd
                   end do
                end do
             end if
          else ! not broyden, use exact jacobian and recalc g

             call this%jac(x_n, F_n, x_scale_, J_a)

             g_a = 0._WP
             
             do j=1,nx
                fafs = F_n(j)*F_scale_(j)**2
                do i=1,nx
                   g_a(i) = g_a(i) + J_a(j,i)*fafs
                end do
             end do
             
          end if

          call this%stop_qq(x, x_n, F_n, ret_status, iter_count, itermax, consec_max, term_code, g_a, f_norm_n, F_scale_, x_scale_)

          if(term_code == 2 .AND. .NOT.(restart)) then
             to_secant = .TRUE.
          else if(term_code > 0) then
             ! done or fatal
             !!print *, 'done or fatal', term_code

             x = x_n
             if(term_code == 1) then
                F_n = F_a
             end if
             f_norm_a = f_norm_n
             F_a = F_n
             no_sol = .FALSE.

          else
             ! progress with algorithm
             restart = .FALSE.
             x = x_n
             f_norm_a = f_norm_n
             F_a = F_n

          end if
          
       end if

    end do
  
    x_out = x

    !!print *, 'iter_count', iter_count
    
    return

    print *, 'never gonna give you up'
    
  end function root_r_
  
  !****

  subroutine jac_r_ (this, x, F_x, x_scale_, J_)

    class(multi_func_t), intent(inout) :: this
    real(WP), intent(inout)           :: x(:)
    real(WP), intent(in)              :: F_x(:)
    real(WP), intent(in)              :: x_scale_(:)
    real(WP), intent(inout)           :: J_(:,:)

    real(WP), allocatable  :: F_j(:)
    
    integer  :: i
    integer  :: j
    integer  :: nx
    integer  :: nxt
    real(WP) :: noise
    real(WP) :: step
    real(WP) :: tempj
    
    nx = SIZE(x)
    
    $ASSERT(nx == this%neq, Wrong size of input array)
    $ASSERT(SIZE(J_, DIM=1)==nx, wrong size of J)
    $ASSERT(SIZE(J_, DIM=2)==nx, wrong size of J)
    
    allocate(F_j(nx))
    
    noise = SQRT(EPSILON(0._WP))

    do j=1,nx

       ! calculate column j of jacobian

       step = noise*MAX(abs(x(j)), 1._WP/x_scale_(j)) * SIGN(1._WP, x(j))
       tempj = x(j)
       x(j) = x(j) + step
       step = x(j) - tempj

       call this%eval_r_(x, F_j)

       do i=1,nx
          J_(i,j) = (F_j(i) - F_x(i))/step
       end do
       x(j) = tempj ! reset
    end do
              
    return
    
  end subroutine jac_r_


  !****

  function gradient_(this, x, f_x, x_scale_) result (del_f)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(inout)           :: x(:)
    real(WP), intent(in)              :: f_x(:)
    real(WP), intent(in)              :: x_scale_(:)

    real(WP), allocatable  :: del_f(:)
    real(WP)               :: f_j
    real(WP)               :: noise

    integer :: nx
    integer :: nxt
    
    nx = SIZE(x)
    $ASSERT(nx == this%neq, Wrong size of input array)

    allocate(del_f(nx))
    
    noise = SQRT(EPSILON(0._WP))

    $ABORT(Need to double check gradient function)
    
    return
    
  end function gradient_
  
  !****

  ! evaluates norm of F_x
  
  function eval_norm_ (this, a, a_scale_) result (f_n)

    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: a(:) 
    real(WP), intent(in), optional    :: a_scale_(:)

    real(WP)             :: f_n
    
    integer :: i
    integer :: nx

    real(WP), allocatable :: a_scale(:)
    
    nx = SIZE(a)
    $ASSERT(nx == this%neq, Wrong size of input array)

    allocate(a_scale(nx))
    
    if(PRESENT(a_scale_)) then
       a_scale = a_scale_
    else
       a_scale = 1._WP
    end if
    
    f_n = 0._WP

    do i=1,nx
       f_n = f_n + (a_scale(i)*a(i))**2
    end do

    f_n = SQRT(f_n)

    return
    
  end function eval_norm_


  ! ****

  subroutine line_search_(this, x_a, x_n, f_norm_a, f_norm_n, F_n, g, p, status, max_taken, F_scale_, x_scale_)
    class(multi_func_t), intent(inout)   :: this
    real(WP), intent(in)                 :: x_a(:)
    real(WP), intent(inout)              :: x_n(:)
    real(WP), intent(in)                 :: f_norm_a  ! from eval_norm
    real(WP), intent(out)                :: f_norm_n
    real(WP), intent(inout)              :: F_n(:)
    real(WP), intent(in)                 :: g(:)    ! 
    real(WP), intent(in)                 :: p(:)    ! Newton step
    integer, intent(out)                 :: status
    logical, intent(out)                 :: max_taken
    real(WP), intent(in)                 :: F_scale_(:)
    real(WP), intent(in)                 :: x_scale_(:)

    real(WP) :: newtlen
    real(WP) :: alpha
    real(WP) :: max_step
    real(WP), allocatable :: p_(:)
    real(WP) :: lambda
    real(WP) :: term_zero
    real(WP) :: term_one
    real(WP) :: term_two
    real(WP) :: a
    real(WP) :: b
    real(WP) :: disc
    real(WP) :: f_n_prev
    real(WP) :: init_slope
    real(WP) :: lambda_prev
    real(WP) :: lambda_temp
    real(WP) :: min_lambda
    real(WP) :: rel_length
    integer  :: nx
    integer  :: i
    integer  :: eval_status
    
    nx = SIZE(x_a)
    $ASSERT(nx == this%neq, Wrong size of input array x)
   
    nx = SIZE(g)
    $ASSERT(nx == this%neq, Wrong size of input array g)

    nx = SIZE(p)
    $ASSERT(nx == this%neq, Wrong size of input array p)

    nx = SIZE(x_n)
    $ASSERT(nx == this%neq, Wrong size of input array p)

    allocate(p_(nx))
    
    alpha = 0.0001_WP

    newtlen = this%eval_norm(p, x_scale_)

    max_step = MAXSTEP_FAC*this%eval_norm(x_a,x_scale_)
    !max_step = MAXSTEP_FAC*MAX(this%eval_norm(x_a,x_scale_), this%eval_norm(x_scale_))
    max_taken = .FALSE.

    p_ = p

    if(newtlen > max_step) then
       p_ = p_*(max_step/newtlen)
       newtlen = max_step
    end if

    init_slope = SUM(g*p_)

    rel_length = 0._WP

    do i=1,nx
       rel_length = MAX(rel_length, ABS(p_(i))/MAX(ABS(x_a(i)), ABS(1._WP/x_scale_(i))))
    end do

    min_lambda = STEP_TOL/rel_length

    ! begin double-checking step
    
    lambda = 1._WP

    status = 2

    do while(status > 1)
      
       x_n = x_a + lambda*p

       call this%eval_r_(x_n, F_n, eval_status)
       f_norm_n = 0.5*this%eval_norm(F_n, F_scale_)**2

       !print *, 'norm check', f_norm_a, f_norm_n

       if(eval_status == EVAL_OK) then
          if(f_norm_n .LE. f_norm_a + alpha*lambda*init_slope) then

             ! Newton step is good
             status = 0
             
             ! check for maximum step
             if(lambda .EQ. 1._WP .AND. newtlen > 0.99 * max_step) then
                max_taken = .TRUE.
             end if
             
          else if(lambda < min_lambda) then
             
             ! lambda too small, fubar
             status = 1
             x_n = x_a
          else
             ! need to backtrack!
             
             if(lambda .EQ. 1) then

                ! first attempt, quad fit
                lambda_temp = -init_slope / (2._WP*(f_norm_n - f_norm_a - init_slope))
             else
                ! other attempts, cubic fit
                
                term_zero = 1._WP/(lambda-lambda_prev)
                term_one = (f_norm_n - f_norm_a - lambda*init_slope)
                term_two = (f_n_prev - f_norm_a - lambda_prev*init_slope)
                
                a = term_zero*(term_one/lambda**2 - term_two/lambda_prev**2)
                b = term_zero*(-lambda_prev*term_one/lambda**2 + lambda*term_two/lambda_prev**2)
                
                if(a .EQ. 0._WP) then
                   lambda_temp = -init_slope/(2._WP*b)
                else
                   disc = b**2 - 3._WP*a*init_slope
                   lambda_temp = (-b+SQRT(disc))/(3._WP*a)
                end if
                
                if(lambda_temp > 0.5_WP * lambda) then
                   lambda_temp = 0.5*lambda
                end if
             end if

             lambda_prev = lambda
             f_n_prev = f_norm_n

             if(lambda_temp .LE. 0.1*lambda) then
                lambda = 0.1 * lambda
             else
                lambda = lambda_temp
             end if

          end if

       else ! eval_status is bad, try backtracking slightly
          !print *, 'here qq'
          lambda_prev = lambda
          f_n_prev = f_norm_n
          lambda = 0.95_WP * lambda
          
       end if

    end do

    return
          
  end subroutine line_search_

  !****

  subroutine hook_search_ (this, x_a, x_n, f_norm_a, f_norm_n, F_n, g, s_n, L, H, F_scale_, x_scale_, iter_count, delta, delta_prev, mu, phi, phi_prime, ret_status, max_taken)

    class(multi_func_t), intent(inout):: this
    real(WP), intent(in)              :: x_a(:)
    real(WP), intent(inout)           :: x_n(:)
    real(WP), intent(in)              :: f_norm_a
    real(WP), intent(out)             :: f_norm_n
    real(WP), intent(inout)           :: F_n(:)
    real(WP), intent(in)              :: g(:) 
    real(WP), intent(inout)           :: s_n(:)
    real(WP), intent(inout)           :: L(:,:)
    real(WP), intent(inout)           :: H(:,:)
    real(WP), intent(inout)           :: F_scale_(:)
    real(WP), intent(inout)           :: x_scale_(:)
    integer, intent(in)               :: iter_count
    real(WP), intent(inout)           :: delta
    real(WP), intent(inout)           :: delta_prev
    real(WP), intent(inout)           :: mu
    real(WP), intent(inout)           :: phi
    real(WP), intent(inout)           :: phi_prime
    integer,  intent(out)             :: ret_status
    logical,  intent(out)             :: max_taken

    logical :: firsthook
    logical :: newt_taken
    real(WP) :: newtlen
    real(WP) :: alpha
    real(WP) :: beta
    real(WP) :: temp
    real(WP) :: max_step
    real(WP) :: phi_prime_init
    real(WP), allocatable :: s(:)
    real(WP), allocatable :: x_n_prev(:)
    real(WP), allocatable :: F_n_prev(:)
    real(WP) :: f_norm_n_prev
    integer :: nx
    integer :: i
    integer :: j

    nx = SIZE(x_a)

    $ASSERT(nx == this%neq, wrong x size)

    allocate(s(nx))
    allocate(x_n_prev(nx))
    allocate(F_n_prev(nx))
    
    ret_status = 4
    firsthook = .TRUE.
    
    max_step = MAXSTEP_FAC*this%eval_norm(x_a,x_scale_)
    newtlen = this%eval_norm(s_n, x_scale_)
    
    if(iter_count == 1 .OR. delta == -1._WP) then
       mu = 0._WP

       if(delta == -1._WP) then

          alpha = this%eval_norm(g, 1._WP/x_scale_)**2
          beta = 0._WP

          do i=1,nx
             temp = 0._WP
             do j=1,nx
                temp = temp + L(j,i)*g(j)/(x_scale_(j)*x_scale_(j))
             end do
             beta = beta + temp*temp
          end do

          delta = alpha*SQRT(alpha)/beta

          if(delta > max_step) then
             delta = max_step
          end if

       end if
    end if

    do while(ret_status > 1)
       call this%hook_step(g, L, H, s_n, x_scale_, newtlen, delta, delta_prev, mu, phi, phi_prime, phi_prime_init, firsthook, s, newt_taken)
       
       delta_prev = delta
       
       call this%trust_region_update(x_a, f_norm_a, g, L, s, F_scale_, x_scale_, newt_taken, H, delta, ret_status, x_n_prev, F_n_prev, f_norm_n_prev, x_n, f_norm_n, F_n, max_taken)
    end do

    return

  end subroutine hook_search_
          
  !****

  subroutine hook_step_ (this, g, L, H, s_n, x_scale_, newtlen, delta, delta_prev, mu, phi, phi_prime, phi_prime_init, firsthook, s, newt_taken)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: g(:)
    real(WP), intent(inout)           :: L(:,:)
    real(WP), intent(inout)           :: H(:,:)
    real(WP), intent(inout)           :: s_n(:)
    real(WP), intent(inout)           :: x_scale_(:)
    real(WP), intent(in)              :: newtlen
    real(WP), intent(inout)           :: delta
    real(WP), intent(in)              :: delta_prev
    real(WP), intent(inout)           :: mu
    real(WP), intent(inout)           :: phi
    real(WP), intent(inout)           :: phi_prime
    real(WP), intent(inout)           :: phi_prime_init
    logical, intent(inout)            :: firsthook
    real(WP), intent(inout)           :: s(:)
    logical,  intent(out)             :: newt_taken

    real(WP), allocatable :: tempvec(:)
    real(WP), allocatable :: tempvec_(:)
    real(WP) :: hi
    real(WP) :: lo
    real(WP) :: mu_low
    real(WP) :: mu_up
    real(WP) :: step_len
    logical  :: done
    
    integer :: nx
    integer :: i

    nx = SIZE(g)

    $ASSERT(nx == this%neq, wrong size of g)

    allocate(tempvec(nx))
    allocate(tempvec_(nx))
    
    hi = 1.5_WP
    lo = 0.75_WP

    if(newtlen .LE. hi*delta) then
       newt_taken = .TRUE.
       s = s_n
       mu = 0._WP
       delta = MIN(delta, newtlen)
    
    else

       newt_taken = .FALSE.

       if(mu > 0._WP) then
          mu = mu - ((phi+delta_prev)/delta) * (((delta_prev - delta)+phi)/phi_prime)
       end if

       phi = newtlen - delta
       
       if(firsthook) then
          firsthook = .FALSE.
          tempvec = x_scale_*x_scale_*s_n
          tempvec_ = tempvec
          call this%lsolve(tempvec_, L, tempvec)

          phi_prime_init = this%eval_norm(tempvec)
          phi_prime_init = -phi_prime_init**2 / newtlen

       end if

       mu_low = -phi/phi_prime_init
       mu_up = this%eval_norm(g, 1._WP/x_scale_)/delta
       done = .FALSE.

       do while(.NOT.(done))
          if((mu < mu_low) .OR. (mu > mu_up)) then
             mu = MAX(SQRT(mu_low*mu_up), 0.001_WP*mu_up)
          end if

          do i=1,nx
             H(i,i) = H(i,i) + mu*x_scale_(i)*x_scale_(i)
          end do

          call this%choldecomp_pd(H, L)
          call this%cholsolve(g, L, s)

          do i=1,nx
             H(i,i) = H(i,i) - mu*x_scale_(i)*x_scale_(i)
          end do

          step_len = this%eval_norm(s, x_scale_)

          phi = step_len - delta

          tempvec_ = s*x_scale_*x_scale_

          call this%lsolve(tempvec_, L, tempvec)

          phi_prime = this%eval_norm(tempvec)
          phi_prime = -phi_prime**2/step_len

          if(((step_len .GE. lo*delta) .AND. (step_len .LE. hi*delta)) .OR. (mu_up - mu_low .LE. 0._WP)) then
             done = .TRUE.
          else
             mu_low = MAX(mu_low, mu - (phi/phi_prime))
             if(phi < 0._WP) then
                mu_up = mu
             end if
             mu = mu - (step_len/delta)*(phi/phi_prime)
          end if
       end do
    end if

    return
    
  end subroutine hook_step_

  !****

  subroutine trust_region_update_ (this, x_a, f_norm_a, g, L, s, F_scale_, x_scale_, newt_taken, H, delta, ret_status, x_n_prev, F_n_prev, f_norm_n_prev, x_n, f_norm_n, F_n, max_taken)

    class(multi_func_t), intent(inout):: this
    real(WP), intent(in)              :: x_a(:)
    real(WP), intent(in)              :: f_norm_a
    real(WP), intent(in)              :: g(:)
    real(WP), intent(in)              :: L(:,:)
    real(WP), intent(in)              :: s(:)
    real(WP), intent(in)              :: F_scale_(:)
    real(WP), intent(in)              :: x_scale_(:)
    logical, intent(in)               :: newt_taken
    real(WP), intent(in)              :: H(:,:)
    real(WP), intent(inout)           :: delta
    integer, intent(inout)            :: ret_status
    real(WP), intent(inout)           :: x_n_prev(:)
    real(WP), intent(inout)           :: F_n_prev(:)
    real(WP), intent(inout)           :: f_norm_n_prev
    real(WP), intent(inout)           :: x_n(:)
    real(WP), intent(inout)           :: f_norm_n
    real(WP), intent(inout)           :: F_n(:)
    logical, intent(out)              :: max_taken

    integer :: nx
    integer :: i
    integer :: j

    real(WP) :: alpha
    real(WP) :: step_len
    real(WP) :: del_f
    real(WP) :: initslope
    real(WP) :: rel_length
    real(WP) :: temp
    real(WP) :: delta_temp
    real(WP) :: del_f_pred
    real(WP) :: max_step
    logical  :: condition
    
    nx = SIZE(x_a)

    $ASSERT(nx == this%neq, wrong size of input array x)

    max_taken = .FALSE.
    max_step = MAXSTEP_FAC*this%eval_norm(x_a,x_scale_)
    
    alpha = 0.0001_WP

    step_len = this%eval_norm(s, x_scale_)

    x_n = x_a + s
    
    call this%eval_r_(x_n, F_n)
    f_norm_n = 0.5*this%eval_norm(F_n, F_scale_)**2

    del_f = f_norm_n - f_norm_a

    initslope = SUM(g*s)

    if (ret_status .NE. 3) then
       f_norm_n_prev = 0._WP
    end if

    if((ret_status == 3) .AND. ((f_norm_n .GE. f_norm_n_prev) .OR. (del_f .GT. alpha*initslope))) then
       ret_status = 0
       x_n = x_n_prev
       f_n = f_n_prev
       delta = delta/2._WP
       F_n = F_n_prev

    elseif (del_f > alpha*initslope) then

       rel_length = 0._WP
       do i=1,nx
          temp = ABS(s(i))/MAX(ABS(x_n(i)), 1._WP/x_scale_(i))
          rel_length = MAX(rel_length, temp)
       end do

       if(rel_length < STEP_TOL) then
          ret_status = 1
          x_n = x_a
       else
          ret_status = 2
          delta_temp = (-initslope*step_len)/(2._WP*(del_f-initslope))

          if(delta_temp < 0.1_WP * delta) then
             delta = 0.1_WP * delta
          else if(delta_temp > 0.5_WP*delta) then
             delta = 0.5_WP*delta
          else
             delta = delta_temp
          end if
       end if
    else
       del_f_pred = initslope

       if(this%global == HOOK_SEARCH) then
          do i=1,nx
             temp = 0.5_WP*H(i,i)*s(i)**2
             do j=i+1,nx
                temp = temp + H(i,j)*s(i)*s(j)
             end do
             del_f_pred = del_f_pred + temp
          end do
       else if(this%global == DOGLEG_SEARCH) then
          do i=1,nx
             temp = 0._WP
             do j=i,nx
                temp = temp + L(j,i)*s(j)
             end do

             del_f_pred = del_f_pred + (temp*temp/2._WP)
          end do
          $ABORT(dogleg not implemented yet)
       else
          $ABORT(what option is this???)
       end if

       ! stupid preprocesser kills the giant if statement condition
       ! had to break it up

       condition = (ret_status .NE. 2) .AND. ((ABS(del_f_pred-del_f) .LE. 0.1_WP*ABS(del_f)) .OR. (del_f .LE. initslope))
       condition = condition .AND. (.NOT.(newt_taken)) .AND. (delta .LE. 0.99_WP*max_step)
       
       if(condition) then

             ! we may be able to take a bigger delta
             ret_status = 3
             x_n_prev = x_n
             f_norm_n_prev = f_norm_n
             F_n_prev = F_n
             delta = MIN(2._WP*delta, max_step)
       else
          ! clear!

          ret_status = 0
          if(step_len > 0.99_WP * max_step) then
             max_taken = .TRUE.
          end if

          if(del_f .GE. 0.1_WP * del_f_pred) then
             delta = delta/2._WP
          else if(del_f .LE. 0.75_WP * del_f_pred) then
             delta = MIN(2._WP*delta, max_step)
          end if

       end if
    end if

    return
    
  end subroutine trust_region_update_
       
  !****

  subroutine nonlinear_solve_(this, restart, F_a, g, M, M2, Jac, H, s_n, F_scale_, x_scale_)
    class(multi_func_t), intent(in)   :: this
    logical, intent(in)               :: restart
    real(WP), intent(in)              :: F_a(:)
    real(WP), intent(in)              :: g(:)
    real(WP), intent(inout)           :: M(:,:)
    real(WP), intent(inout)           :: M2(:)
    real(WP), intent(inout)           :: Jac(:,:)
    real(WP), intent(inout)           :: H(:,:)
    real(WP), intent(inout)           :: s_n(:)
    real(WP), intent(in)              :: F_scale_(:)
    real(WP), intent(in)              :: x_scale_(:)
    
    real(WP), allocatable :: M1(:)

    real(WP) :: est
    real(WP) :: fsf
    real(WP) :: Hnorm
    real(WP) :: temp
    
    logical :: singular
    
    integer :: nx
    integer :: i
    integer :: j
    integer :: k

    nx = SIZE(F_a)

    $ASSERT(nx == this%neq, Wrong size of input array F_a)

    $ASSERT(SIZE(M2) == nx, wrong size of array M2)
    $ASSERT(SIZE(M, DIM=2) == nx, wrong size of array M)
    $ASSERT(SIZE(M, DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(Jac, DIM=2) == nx, wrong size of array Jac)
    $ASSERT(SIZE(Jac, DIM=1) == nx, wrong size of array Jac)

    allocate(M1(nx))

    if(restart) then
       do i=1,nx
          do j=1,nx
             M(i,j) = F_scale_(i)*Jac(i,j)
          end do
       end do

       call this%QR_decomp(M, M1, M2, singular)

       call this%Q_form(M,M1,Jac)

       do i=1,nx
          M(i,i) = M2(i)
       end do

    else

       singular = .FALSE.

       do i=1,nx
          if(M(i,i) == 0._WP) then
             singular = .TRUE.
          end if
       end do
    end if
    

    if(.NOT.(singular)) then
       do j=1,nx
          do i=1,j-1
             M(i,j) = M(i,j)/x_scale_(j)
          end do
          M2(j) = M2(j)/x_scale_(j)
       end do
       
       call this%condest(M,M2,est)
       
       do j=1,nx
          do i=1,j-1
             M(i,j) = M(i,j)*x_scale_(j)
          end do
          M2(j) = M2(j)*x_scale_(j)
       end do
          
    else !if singular
          
       est = 0._WP
       
    end if
    
    if(singular .OR. est > 1._WP/sqrt(EPSILON(0._WP))) then

       do i=1,nx
          do j=1,nx
             H(i,j) = 0._WP
             do k=1,i
                H(i,j) = H(i,j) + M(k,i)*M(k,j)
             end do
          end do
       end do
       
       Hnorm = 0._WP
       do j=1,nx
          Hnorm = Hnorm + ABS(H(1,j))/x_scale_(j)
       end do
       
       Hnorm = Hnorm/x_scale_(1)
       
       do i=2,nx
          temp = 0._WP
          do j=1,i
             temp = temp + ABS(H(j,i))/x_scale_(j)
          end do
          do j=(i+1),nx
             temp = temp + ABS(H(i,j))/x_scale_(j)
          end do
          
          temp = temp/x_scale_(i)
          
          Hnorm = MAX(temp, Hnorm)
       end do
       
       do i=1,nx
          H(i,i) = H(i,i) + SQRT(nx*EPSILON(0._WP)) * Hnorm * x_scale_(i)**2
       end do
       
       call this%choldecomp_pd(H, M)
       call this%cholsolve(g,M,s_n)
       
    else

       s_n = 0._WP
       do j=1,nx
          fsf = F_scale_(j)*F_a(j)
          do i=1,nx
             s_n(i) = s_n(i) - Jac(i,j)*fsf
          end do
       end do
       
       call this%rsolve(M, M2, s_n)
       
    end if
 

    
    return  
       
  end subroutine nonlinear_solve_

  !****

  subroutine nonlinear_solve_unfac_(this, F_a, g, Jac, M, H, s_n, F_scale_, x_scale_)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: F_a(:)
    real(WP), intent(in)              :: g(:)
    real(WP), intent(in)              :: Jac(:,:)
    real(WP), intent(inout)           :: M(:,:)
    real(WP), intent(inout)           :: H(:,:)
    real(WP), intent(inout)           :: s_n(:)
    real(WP), intent(in)              :: F_scale_(:)
    real(WP), intent(in)              :: x_scale_(:)

    integer :: nx
    integer :: i
    integer :: j
    integer :: k
    logical :: singular

    real(WP) :: est
    real(WP) :: Hnorm
    real(WP) :: temp
    real(WP), allocatable :: M1(:)
    real(WP), allocatable :: M2(:)
    
    ! reset M, H, s_n since they are really "output" variables

    M = 0._WP
    H = 0._WP
    s_n = 0._WP
    
    nx = SIZE(F_a)

    $ASSERT(nx == this%neq, Wrong size of input array F_a)

    $ASSERT(SIZE(M, DIM=2) == nx, wrong size of array M)
    $ASSERT(SIZE(M, DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(Jac, DIM=2) == nx, wrong size of array Jac)
    $ASSERT(SIZE(Jac, DIM=1) == nx, wrong size of array Jac)
    $ASSERT(SIZE(H, DIM=2) == nx, wrong size of array H)
    $ASSERT(SIZE(H, DIM=1) == nx, wrong size of array H)
    $ASSERT(SIZE(s_n) == nx, wrong size of array s_n)

    allocate(M1(nx))
    allocate(M2(nx))
    
    do i=1,nx
       do j=1,nx
          M(i,j) = F_scale_(i)*Jac(i,j)
       end do
    end do

    call this%QR_decomp(M, M1,M2,singular)

    if(.NOT.(singular)) then

       do j=1,nx
          do i=1,j-1
             M(i,j) = M(i,j)/x_scale_(j)
          end do
          M2(j) = M2(j)/x_scale_(j)
       end do

       call this%condest(M,M2,est)

    else
       est = 0._WP
    end if

    if(singular .OR. est > 1._WP/SQRT(EPSILON(0._WP))) then

       do i=1,nx
          do j=i,nx
             H(i,j) = 0._WP
             do k=1,nx
                H(i,j) = H(i,j) + Jac(k,i)*Jac(k,j)*F_scale_(k)**2
             end do
          end do
       end do

       Hnorm = 0._WP
       do j=1,nx
          Hnorm = Hnorm + ABS(H(1,j))/x_scale_(j)
       end do

       Hnorm = Hnorm/x_scale_(1)
       
       do i=2,nx
          temp = 0._WP
          do j=1,i
             temp = temp + ABS(H(j,i))/x_scale_(j)
          end do
          do j=(i+1),nx
             temp = temp + ABS(H(i,j))/x_scale_(j)
          end do
          
          temp = temp/x_scale_(i)
          
          Hnorm = MAX(temp, Hnorm)
       end do

       do i=1,nx
          H(i,i) = H(i,i) + SQRT(nx*EPSILON(0._WP)) * Hnorm * x_scale_(i)**2
       end do

       call this%choldecomp_pd(H, M)
       call this%cholsolve(g,M,s_n)

    else

       do j=1,nx
          do i=1,j-1
             M(i,j) = M(i,j)*x_scale_(j)
          end do
          M2(j) = M2(j)*x_scale_(j)
       end do

       s_n = -F_a*F_scale_

       call this%QR_solve(M,M1,M2,s_n)

       if(this%global .NE. LINE_SEARCH) then
          do i=1,nx
             M(i,i) = M2(i)
             do j=1,i-1
                M(i,j) = M(j,i)
             end do
          end do
       end if

       if(this%global == HOOK_SEARCH) then
          do i=1,nx
             H(i,i) = 0._WP
             do k=1,i
                H(i,i) = H(i,i) + M(i,k)**2
             end do
             do j=i+1,nx
                H(i,j) = 0._WP
                do k=1,i
                   H(i,j) = H(i,j) + M(i,k)*M(j,k)
                end do
             end do

          end do
       end if
    end if
    
    return
    
  end subroutine nonlinear_solve_unfac_
    
  !****

  subroutine QR_solve_(this, M, M1, M2, b)
    class(multi_func_t), intent(in)    :: this
    real(WP), intent(in)               :: M(:,:)
    real(WP), intent(in)               :: M1(:)
    real(WP), intent(in)               :: M2(:)
    real(WP), intent(inout)            :: b(:)

    integer :: i
    integer :: j
    integer :: nx
    real(WP) :: tau
    
    nx = SIZE(M1)
    $ASSERT(nx == this%neq, Wrong size of input array M1)
    $ASSERT(SIZE(M, DIM=2) == nx, wrong size of array M)
    $ASSERT(SIZE(M, DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(M2) == nx, wrong size of array M2)
    $ASSERT(SIZE(b) == nx, wrong size of array b)

    do j=1,nx-1
       tau = 0._WP
       do i=j,nx
          tau = tau + M(i,j)*b(i)
       end do
       tau = tau/M1(j)

       do i=j,nx
          b(i) = b(i) - tau*M(i,j)
       end do

    end do

    call this%rsolve(M,M2,b)

    return
    
  end subroutine QR_solve_
    

  !****

  subroutine QR_decomp_(this, M, M1, M2, singular)
    class(multi_func_t), intent(in)    :: this
    real(WP), intent(inout)            :: M(:,:)
    real(WP), intent(inout)            :: M1(:)
    real(WP), intent(inout)            :: M2(:)
    logical, intent(out)               :: singular

    integer :: i
    integer :: j
    integer :: k
    integer :: nx

    real(WP) :: eta
    real(WP) :: sigma
    real(WP) :: tau
    
    singular = .FALSE.

    nx = this%neq

    $ASSERT(SIZE(M, DIM=1) == nx, wrong size of M)
    $ASSERT(SIZE(M, DIM=2) == nx, wrong size of M)
    $ASSERT(SIZE(M1) == nx, wrong size of M1)
    $ASSERT(SIZE(M2) == nx, wrong size of M2)

    do k=1,(nx-1)
       eta = 0._WP

       do i=k,nx
          eta = MAX(eta, ABS(M(i,k)))
       end do

       if(eta == 0._WP) then
          singular = .TRUE.
          M1(k) = 0._WP
          M2(k) = 0._WP

       else
          sigma = 0._WP
          
          do i =k,nx
             M(i,k) = M(i,k) / eta
             sigma = sigma + M(i,k)**2
          end do

          sigma = SQRT(sigma) * SIGN(1._WP, M(k,k))

          M(k,k) = M(k,k) + sigma
          M1(k) = sigma*M(k,k)
          M2(k) = -eta*sigma

          do j=(k+1),nx
             tau = 0._WP
             do i=k,nx
                tau = tau+M(i,k)*M(i,j)
             end do
             tau = tau/M1(k)

             do i=k,nx
                M(i,j) = M(i,j) - tau*M(i,k)
             end do
          end do
       end if
    end do
    
    if(M(nx,nx) == 0._WP) then
       singular = .TRUE.
    end if
    
    M2(nx) = M(nx,nx)

    return       
   
  end subroutine QR_decomp_

  !****

  subroutine QR_update_(this, u, v, Z, M)

    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: u(:)
    real(WP), intent(in)              :: v(:)
    real(WP), intent(inout)           :: Z(:,:)
    real(WP), intent(inout)           :: M(:,:)

    integer :: nx
    integer :: i
    integer :: j
    integer :: k

    real(WP), allocatable :: u_(:)

    nx = this%neq

    $ASSERT(SIZE(Z,DIM=1) == nx, wrong size of array Z)
    $ASSERT(SIZE(Z,DIM=2) == nx, wrong size of array Z)
    $ASSERT(SIZE(M,DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(M,DIM=2) == nx, wrong size of array M)
    
    allocate(u_(nx))

    u_ = u
    
    do i=2,nx
       M(i,i-1) = 0._WP
    end do

    k = nx

    do while((u_(k) .EQ. 0._WP) .AND. (k > 1))
       k = k-1
    end do

    do i=(k-1),1,-1
       call this%jac_rotate(i, u_(i), -u_(i+1), Z, M)

       if(u_(i) .EQ. 0._WP) then
          u_(i) = ABS(u_(i+1))
       else
          u_(i) = SQRT(u_(i)**2 + u_(i+1)**2)
       end if

    end do

    do j=1,nx
       M(1,j) = M(1,j) + u_(1)*v(j)
    end do

    do i=1,k-1
       call this%jac_rotate(i, M(i,i), -M(i+1,i), Z,M)
    end do

    return
    
  end subroutine QR_update_
  
    
  !****

  subroutine Q_form_(this, M, M1, Z)

    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: M(:,:)
    real(WP), intent(in)              :: M1(:)
    real(WP), intent(inout)           :: Z(:,:)

    integer :: nx
    integer :: i
    integer :: j
    integer :: k

    real(WP) :: tau
    
    nx = this%neq
    $ASSERT(SIZE(Z, DIM=1) == nx, wrong size of Z)
    $ASSERT(SIZE(Z, DIM=2) == nx, wrong size of Z)

    Z = 0._WP
    do i=1,nx
       Z(i,i) = 1._WP
    end do

    do k=1,nx-1
       if(M1(k) .NE. 0._WP) then
          do j=1,nx
             tau = 0._WP
             do i=k,nx
                tau = tau + M(i,k)*Z(i,j)
             end do
             tau = tau/M1(k)

             do i=k,nx
                Z(i,j) = Z(i,j) - tau*M(i,k)
             end do
          end do
       end if
    end do
    
    return

  end subroutine Q_form_
  !****

  subroutine condest_(this,M,M2,est)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: M(:,:)
    real(WP), intent(in)              :: M2(:)
    real(WP), intent(out)             :: est

    real(WP), allocatable :: p(:)
    real(WP), allocatable :: pm(:)
    real(WP), allocatable :: x(:)

    real(WP) :: temp
    real(WP) :: tempm
    real(WP) :: xp
    real(WP) :: xm
    real(WP) :: xnorm
    
    integer :: nx
    integer :: i
    integer :: j
    integer :: k

    nx = this%neq

    allocate(p(nx))
    allocate(pm(nx))
    allocate(x(nx))

    est = ABS(M2(1))

    do j=2,nx
       temp = ABS(M2(j))
       do i=1,j-1
          temp = temp+ABS(M(i,j))
       end do
       est = MAX(temp, est)
    end do

    x(1)  = 1._WP / M2(1)

    do i =2,nx
       p(i) = M(1,i)*x(1)
    end do
    
    
    do j=2,nx
       xp = (1._WP - p(j)) / M2(j)
       xm = (-1._WP - p(j)) / M2(j)
       temp = ABS(xp)
       tempm = ABS(xm)

       do i=(j+1),nx
          pm(i) = p(i) +M(j,i)*xm
          tempm=tempm+(ABS(pm(i))/ABS(M2(i)))
          p(i) = p(i) + M(j,i)*xp
          temp = temp+(ABS(p(i))/ABS(M2(i)))
       end do

       if(temp > tempm) then
          x(j) = xp
       else
          x(j) = xm
          do i=(j+1),nx
             p(i) = pm(i)
          end do
       end if
    end do

    xnorm = SUM(ABS(x))

    est = est / xnorm

    call this%rsolve(M,M2,x)

    xnorm = SUM(ABS(x))

    est = est*xnorm

    return
    
  end subroutine condest_
  
  !****

  subroutine rsolve_(this, M, M2, b)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: M(:,:)
    real(WP), intent(in)              :: M2(:)
    real(WP), intent(inout)           :: b(:)

    integer :: nx
    integer :: i
    integer :: j

    real(WP) :: temp
    
    nx = this%neq

    $ASSERT(SIZE(b) == nx, wrong size of array b)
    
    b(nx) = b(nx)/M2(nx)

    do i=(nx-1),1,-1
       temp = 0._WP
       do j=(i+1),nx
          temp = temp + M(i,j)*b(j)
       end do

       b(i) = (b(i) - temp)/M2(i)

    end do

    return       
    
  end subroutine rsolve_
    
  !****

  subroutine cholsolve_(this, g, L, s)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: g(:)
    real(WP), intent(in)              :: L(:,:)
    real(WP), intent(inout)           :: s(:)

    integer :: nx
    integer :: i
    integer :: j

    real(WP), allocatable :: s_temp(:)
    
    nx = this%neq

    $ASSERT(SIZE(s) == nx, wrong size of array s)
    
    allocate(s_temp(nx))
    
    call this%lsolve(g,L,s)
    call this%ltsolve(s,L,s_temp)

    $ASSERT(SIZE(s) == SIZE(s_temp), something weird broke in lsolve/ltsolve)
    
    s = -s_temp

    return

  end subroutine cholsolve_
    
  !****

  ! H is known to be positive definitive

  subroutine choldecomp_pd_(this, H, L)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: H(:,:)
    real(WP), intent(inout)           :: L(:,:)

    integer :: nx
    integer :: i
    integer :: j
    integer :: k
    real(WP) :: maxoffl
    real(WP) :: minl
    real(WP) :: minl2
    real(WP) :: minljj
    real(WP) :: temp
    real(WP)  :: maxadd    
    
    nx = this%neq
    
    $ASSERT(SIZE(H,DIM=1) == nx, wrong size of array H)
    $ASSERT(SIZE(H,DIM=2) == nx, wrong size of array H)
    $ASSERT(SIZE(L,DIM=1) == nx, wrong size of array L)
    $ASSERT(SIZE(L,DIM=2) == nx, wrong size of array L)

    minl = 0._WP
    maxadd = 0._WP

    maxoffl = 0._WP
    do i=1,nx
       maxoffl = MAX(ABS(H(i,i)), maxoffl)
    end do
    maxoffl = SQRT(maxoffl)

    minl2 = SQRT(EPSILON(0._WP)) * maxoffl

    do j=1,nx
       temp = 0._WP
       do i=1,j-1
          temp = temp + L(j,i)**2
       end do
              
       L(j,j) = H(j,j) - temp

       minljj = 0._WP

       do i=(j+1),nx
          temp = 0._WP
          do k=1,j-1
             temp = temp + L(i,k)*L(j,k)
          end do

          L(i,j) = H(j,i) - temp
          minljj = MAX(ABS(L(i,j)), minljj)
       end do

       minljj = MAX(minl, minljj/maxoffl)

       if(L(j,j) > minljj**2) then
          L(j,j) = sqrt(L(j,j))
       else
          if(minljj < minl2) then
             minljj = minl2
          end if
          maxadd = MAX(maxadd, minljj**2 - L(j,j))
          L(j,j) = minljj
       end if

       do i=(j+1),nx
          L(i,j) = L(i,j)/L(j,j)
       end do
    end do

    return
    
  end subroutine choldecomp_pd_
             
       
  !****

  subroutine lsolve_(this, b, L, y)
    
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: b(:)
    real(WP), intent(in)              :: L(:,:)
    real(WP), intent(inout)           :: y(:)

    integer :: nx
    integer :: i
    integer :: j
    real(WP) :: temp
    
    nx = this%neq

    $ASSERT(SIZE(y) == nx, wrong size of y array)
    
    y(1) = b(1)/L(1,1)

    do i=2,nx
       temp = 0._WP
       do j=1,i-1
          temp = temp + L(i,j)*y(j)
       end do
       y(i) = (b(i)-temp)/L(i,i)
    end do
   
    return
    
  end subroutine lsolve_
  !****

  subroutine ltsolve_(this, y, L,x)
    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: y(:)
    real(WP), intent(in)              :: L(:,:)
    real(WP), intent(inout)             :: x(:)

    integer :: nx
    integer :: i
    integer :: j
    real(WP) :: temp
    
    nx = this%neq

    $ASSERT(SIZE(x) == nx, wrong size of x array)

    x(nx) = y(nx)/L(nx,nx)

    do i=(nx-1),1,-1
       temp = 0._WP
       do j=(i+1),nx
          temp = temp + L(j,i)*x(j)
       end do

       x(i) = (y(i) - temp)/L(i,i)
    end do

    return

  end subroutine ltsolve_

  !****

  subroutine broyden_fac_(this, x_a, x_n, F_a, F_n, Z, M, M2, F_scale_, x_scale_)

    class(multi_func_t), intent(in)   :: this
    real(WP), intent(in)              :: x_a(:)
    real(WP), intent(in)              :: x_n(:)
    real(WP), intent(in)              :: F_a(:)
    real(WP), intent(in)              :: F_n(:)
    real(WP), intent(inout)           :: Z(:,:)
    real(WP), intent(inout)           :: M(:,:)
    real(WP), intent(inout)           :: M2(:)
    real(WP), intent(in)              :: F_scale_(:)
    real(WP), intent(in)              :: x_scale_(:)

    real(WP), allocatable :: s(:)
    real(WP), allocatable :: t(:)
    real(WP), allocatable :: w(:)

    real(WP) :: temp
    real(WP) :: denom
    logical :: skip_update
    integer :: nx
    integer :: i
    integer :: j

    nx = this%neq

    $ASSERT(SIZE(Z,DIM=1) == nx, wrong size of array Z)
    $ASSERT(SIZE(Z,DIM=2) == nx, wrong size of array Z)
    $ASSERT(SIZE(M,DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(M,DIM=2) == nx, wrong size of array M)
    $ASSERT(SIZE(M2) == nx, wrong size of array M2)
    
    allocate(s(nx))
    allocate(t(nx))
    allocate(w(nx))

    do i=1,nx
       M(i,i) = M2(i)
    end do

    s = x_n - x_a
    
    skip_update = .TRUE.

    t = 0._WP
    do i=1,nx
       do j=1,nx
          t(i) = t(i) + M(i,j)*s(j)
       end do
    end do

    do i=1,nx
       temp = 0._WP
       do j=1,nx
          temp = temp + Z(j,i)*t(j)
       end do

       w(i) = F_scale_(i)*(F_n(i) - F_a(i)) - temp

       if(w(i) > EPSILON(0._WP)*F_scale_(i)*(ABS(F_n(i)) - ABS(F_a(i)))) then
          skip_update = .FALSE.
       else
          w(i) = 0._WP
       end if
    end do
    
    if(.NOT.(skip_update)) then
       do i=1,nx
          t(i) = 0._WP
          do j=1,nx
             t(i) = t(i) + Z(i,j)*w(j)
          end do
       end do
       
       denom  = 0._WP
       do i=1,nx
          denom = denom + (x_scale_(i)*s(i))**2
       end do

       denom = SQRT(0.5_WP*denom)
       !print *, 'denom', denom, s
       s = x_scale_**2 * s / denom

       call this%QR_update(t,s,Z,M)
       do i=1,nx
          M2(i) = M(i,i)
       end do
    end if

    return
    
  end subroutine broyden_fac_

  !****

  subroutine jac_rotate_(this, i, a, b, Z, M)
    class(multi_func_t), intent(in)    :: this
    integer, intent(in)                :: i
    real(WP), intent(in)               :: a
    real(WP), intent(in)               :: b
    real(WP), intent(inout)            :: Z(:,:)
    real(WP), intent(inout)            :: M(:,:)

    real(WP) :: c
    real(WP) :: s
    real(WP) :: den
    real(WP) :: y
    real(WP) :: w

    integer :: nx
    integer :: j
    
    nx = this%neq


    $ASSERT(SIZE(Z,DIM=1) == nx, wrong size of array Z)
    $ASSERT(SIZE(Z,DIM=2) == nx, wrong size of array Z)
    $ASSERT(SIZE(M,DIM=1) == nx, wrong size of array M)
    $ASSERT(SIZE(M,DIM=2) == nx, wrong size of array M)
    
    if(a == 0._WP) then
       c = 0._WP
       s = SIGN(1._WP, b)
    else
       den = SQRT(a**2 + b**2)
       c = a/den
       s = b/den
    end if

    do j=1,nx
       y = M(i,j)
       w = M(i+1,j)
       M(i,j) = c*y - s*w
       M(i+1,j) = s*y + c*w

       y = Z(i,j)
       w = Z(i+1,j)
       Z(i,j) = c*y - s*w
       Z(i+1,j) = s*y + c*w
    end do

    return
    
  end subroutine jac_rotate_

  !****

  subroutine stop_qq_(this, x_a, x_n, F_n, ret_status, itncount, itnlimit, consec_max, term_code, g_a, f_norm_n, F_scale_, x_scale_)
    class(multi_func_t), intent(in) :: this
    real(WP), intent(in)     :: x_a(:)
    real(WP), intent(in)     :: x_n(:)
    real(WP), intent(in)     :: F_n(:)
    integer, intent(in)      :: ret_status
    integer, intent(in)      :: itncount
    integer, intent(in)      :: itnlimit
    integer, intent(in)      :: consec_max
    integer, intent(out)     :: term_code
    real(WP), intent(in)     :: g_a(:)
    real(WP), intent(in)     :: f_norm_n
    real(WP), intent(in)     :: F_scale_(:)
    real(WP), intent(in)     :: x_scale_(:)

    integer :: nx
    integer :: i
    
    real(WP) :: max_scaled_func
    real(WP) :: max_scaled_step
    real(WP) :: temp
    real(WP) :: local_min_F
    
    nx = SIZE(x_a)

    term_code = 0

    !!print *, 'iter', itncount, F_n

    max_scaled_func = MAXVAL(F_scale_*ABS(F_n))

    if(max_scaled_func < TOL) then
          term_code = ALL_OK
    else if(ret_status == 1) then
       term_code = NORM_DIDNT_DECREASE
    else if(itncount > itnlimit) then
       term_code = MAX_ITER_EXCEEDED
    else if(consec_max > 5) then
       term_code = QUITE_LARGE_STEPS
    else if(.NOT.(this%broyden)) then
       local_min_F = 0._WP
       do i=1,nx
          temp = ABS(g_a(i))*MAX(ABS(x_n(i)), 1._WP/x_scale_(i))/ MAX(f_norm_n, REAL(nx)/2._WP)
          local_min_F = MAX(local_min_F, temp)
       end do

       if(local_min_F .LE. NORM_TOL) then
          term_code = MAY_NEED_DIFFERENT_INITIAL
       end if
       
    else
       max_scaled_step = 0._WP
       do i=1,nx
          temp = ABS(x_n(i) - x_a(i))/MAX(ABS(x_n(i)), 1._WP/x_scale_(i))
          max_scaled_step = MAX(max_scaled_step, temp)
       end do

       if(max_scaled_step < STEP_TOL) then
          term_code = TOO_SLOW
       end if
    end if

    return
    
  end subroutine stop_qq_
  
end module core_multi_func
