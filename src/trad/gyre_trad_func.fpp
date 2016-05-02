! Program  : gyre_trad_func
! Purpose  : traditional approximation functions
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

module gyre_trad_func

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel

  use astro_hough

  use gyre_cheby

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: trad_func_t
     private
     type(cheby_t)   :: cb_neg
     type(cheby_t)   :: cb_pos
     type(cheby_t)   :: cb_ctr
     real(WP)        :: nu_t
     integer, public :: m
     integer, public :: k
   contains
     private
     procedure       :: lambda_r_
     procedure       :: lambda_c_
     generic, public :: lambda => lambda_r_, lambda_c_
  end type trad_func_t

  ! Interfaces

  interface trad_func_t
     module procedure trad_func_t_tol_
  end interface trad_func_t

  $if ($HDF5)
  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write
  $endif

  ! Access specifiers

  private

  public :: trad_func_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif

  ! Procedures

contains

  function trad_func_t_tol_ (m, k, lambda_tol, cheby_tol, cheby_n) result (tf)

    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: lambda_tol
    real(WP), intent(in) :: cheby_tol
    integer, intent(in)  :: cheby_n
    type(trad_func_t)    :: tf

    integer :: l

    ! Construct the trad_func_t with the specified tolerances and
    ! Chebyshev fit order

    tf%m = m
    tf%k = k

    if (k >= 0) then

       ! Gravity waves

       tf%nu_t = 1._WP

       tf%cb_neg = cheby_t(-1._WP, 0._WP, cheby_n, f_grav_o_)
       tf%cb_pos = cheby_t(0._WP , 1._WP, cheby_n, f_grav_o_)
       tf%cb_ctr = cheby_t(-1._WP, 1._WP, cheby_n, f_grav_i_)

       call tf%cb_neg%truncate(cheby_tol)
       call tf%cb_pos%truncate(cheby_tol)
       call tf%cb_ctr%truncate(cheby_tol)

    else

       $ASSERT(m /= 0,Invalid m for Rossby waves)

       ! Rossby waves

       l = ABS(m) + ABS(k) - 1
       tf%nu_t = REAL(l*(l+1), WP)/REAL(m)

       tf%cb_neg = cheby_t(-1._WP, 0._WP, cheby_n, f_ross_)

       call tf%cb_neg%truncate(cheby_tol)

    endif
       
    ! Finish

    return

  contains

    function f_grav_i_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the gravity-wave eigenvalue function in the inner
      ! region (|nu| < nu_t)

      if (m == 0._WP .AND. k == 0._WP) then
         
         f = 1._WP

      else

         associate (nu => tf%nu_t*x)
           f = lambda(nu, m, k, lambda_tol)/lambda_norm_grav_i_r_(nu, m, k)
         end associate

      endif

      ! Finish

      return

    end function f_grav_i_

    function f_grav_o_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the gravity-wave eigenvalue function in the outer
      ! region (|nu| > nu_t)

      if (m == 0._WP .AND. k == 0._WP) then

         f = 1._WP

      else

         if (x == 0._WP) then

            f = 1._WP

         else

            associate (nu => tf%nu_t/x)
              f = lambda(nu, m, k, lambda_tol)/lambda_norm_grav_o_r_(nu, m, k)
            end associate

         endif

      endif

      ! Finish

      return

    end function f_grav_o_

    function f_ross_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the Rossby-wave eigenvaule function

      if (x == 0._WP) then

         f = 1._WP

      elseif (x == -1._WP) then

         f = 0._WP

      else

         associate (nu => tf%nu_t/x)
           f = lambda(nu, m, k, lambda_tol)/lambda_norm_ross_r_(nu, m, k)
         end associate

      endif

    end function f_ross_

  end function trad_func_t_tol_

  !****

  $if ($HDF5)

  subroutine read_ (hg, tf)

    type(hgroup_t), intent(inout)  :: hg
    type(trad_func_t), intent(out) :: tf

    type(hgroup_t) :: hg_comp

    ! Read the trad_func_t

    call read_attr(hg, 'm', tf%m)
    call read_attr(hg, 'k', tf%k)

    call read_attr(hg, 'nu_t', tf%nu_t)

    hg_comp = hgroup_t(hg, 'cb_neg')
    call read(hg_comp, tf%cb_neg)
    call hg_comp%final()

    if (tf%k >= 0) then

       hg_comp = hgroup_t(hg, 'cb_pos')
       call read(hg_comp, tf%cb_pos)
       call hg_comp%final()
       
       hg_comp = hgroup_t(hg, 'cb_ctr')
       call read(hg_comp, tf%cb_ctr)
       call hg_comp%final()

    end if

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, tf)

    type(hgroup_t), intent(inout) :: hg
    type(trad_func_t), intent(in) :: tf

    type(hgroup_t) :: hg_comp

    ! Write the trad_func_t

    call write_attr(hg, 'm', tf%m)
    call write_attr(hg, 'k', tf%k)

    call write_attr(hg, 'nu_t', tf%nu_t)

    hg_comp = hgroup_t(hg, 'cb_neg')
    call write(hg_comp, tf%cb_neg)
    call hg_comp%final()

    if (tf%k >= 0) then

       hg_comp = hgroup_t(hg, 'cb_pos')
       call write(hg_comp, tf%cb_pos)
       call hg_comp%final()

       hg_comp = hgroup_t(hg, 'cb_ctr')
       call write(hg_comp, tf%cb_ctr)
       call hg_comp%final()

    endif

    ! Finish

    return

  end subroutine write_

  $endif

  !****

  function lambda_r_ (this, nu) result (lambda)

    class(trad_func_t), intent(in) :: this
    real(WP), intent(in)           :: nu
    real(WP)                       :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (real)

    if (this%k >= 0) then

       ! Gravity waves

       if (nu <= -this%nu_t) then
          lambda = this%cb_neg%eval(this%nu_t/nu)*lambda_norm_grav_o_r_(nu, this%m, this%k)
       elseif (nu >= this%nu_t) then
          lambda = this%cb_pos%eval(this%nu_t/nu)*lambda_norm_grav_o_r_(nu, this%m, this%k)
       else
          lambda = this%cb_ctr%eval(nu/this%nu_t)*lambda_norm_grav_i_r_(nu, this%m, this%k)
       endif

    else

       ! Rossby waves

       if (nu <= -this%nu_t) then
          lambda = this%cb_neg%eval(this%nu_t/nu)*lambda_norm_ross_r_(nu, this%m, this%k)
       else
          $ABORT(Invalid nu for Rossby waves)
       end if

    end if

    ! Finish

    return

  end function lambda_r_

  !****

  function lambda_c_ (this, nu) result (lambda)

    class(trad_func_t), intent(in) :: this
    complex(WP), intent(in)        :: nu
    complex(WP)                    :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (complex)

    if (this%k >= 0) then

       ! Gravity waves

       if (REAL(nu) < -this%nu_t) then
          lambda = this%cb_neg%eval(this%nu_t/nu)*lambda_norm_grav_o_c_(nu, this%m, this%k)
       elseif (REAL(nu) > this%nu_t) then
          lambda = this%cb_pos%eval(this%nu_t/nu)*lambda_norm_grav_o_c_(nu, this%m, this%k)
       else
          lambda = this%cb_ctr%eval(nu/this%nu_t)*lambda_norm_grav_i_c_(nu, this%m, this%k)
       endif

    else

       ! Rossby waves

       if (REAL(nu) <= -this%nu_t) then
          lambda = this%cb_neg%eval(this%nu_t/nu)*lambda_norm_ross_c_(nu, this%m, this%k)
       else
          $ABORT(Invalid nu for Rossby waves)
       end if

    end if

    ! Finish

    return

  end function lambda_c_

  !****

  $define $LAMBDA_NORM $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function lambda_norm_grav_i_${SUFFIX}_ (nu, m, k) result (lambda_norm)

    $TYPE(WP), intent(in) :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    $TYPE(WP)             :: lambda_norm

    ! Evaluate the gravity-wave eigenvalue normalization function in
    ! the inner region (|nu| < nu_t)

    associate (l => ABS(m) + k)
      lambda_norm = l*(l+1)
    end associate

    ! Finish

    return
    
  end function lambda_norm_grav_i_${SUFFIX}_
  
  !****

  function lambda_norm_grav_o_${SUFFIX}_ (nu, m, k) result (lambda_norm)

    $TYPE(WP), intent(in) :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    $TYPE(WP)             :: lambda_norm

    ! Evaluate the gravity-wave eigenvalue normalization function in
    ! the outer region (|nu| > nu_t)

    lambda_norm = lambda_asymp_${SUFFIX}_(nu, m, k)

    ! Finish

    return
    
  end function lambda_norm_grav_o_${SUFFIX}_
  
  !****

  function lambda_norm_ross_${SUFFIX}_ (nu, m, k) result (lambda_norm)

    $TYPE(WP), intent(in) :: nu
    integer, intent(in)   :: m
    integer, intent(in)   :: k
    $TYPE(WP)             :: lambda_norm

    ! Evaluate the Rossby-wave eigenvalue normalization function

    lambda_norm = lambda_asymp_${SUFFIX}_(nu, m, k)

    ! Finish

    return
    
  end function lambda_norm_ross_${SUFFIX}_

  $endsub

  $LAMBDA_NORM(r,real)
  $LAMBDA_NORM(c,complex)
  
  !****

  function lambda_asymp_r_ (nu, m, k) result (lambda)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda

    integer :: s
    
    ! Evaluate the (m, k) Hough eigenvalue lambda using the asymptotic
    ! expressions by [Tow2003a]

    if (m*nu >= 0._WP) then

       $ASSERT(k >= 0,Invalid k)

       if (k > 0) then

          s = k - 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))


       else

          lambda = m**2*(2._WP*m*nu)/(2._WP*m*nu - 1._WP)

       endif

    else

       if (k >= -1) then

          s = k + 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))

       else

          s = -k -1

          lambda = (m*nu - m**2)**2/(nu**2*(2*s+1)**2)

       endif

    endif

    ! Finish

    return

  end function lambda_asymp_r_

  !****

  function lambda_asymp_c_ (nu, m, k) result (lambda)

    complex(WP), intent(in) :: nu
    integer, intent(in)     :: m
    integer, intent(in)     :: k
    complex(WP)             :: lambda

    integer :: s
    
    ! Evaluate the (m, k) Hough eigenvalue lambda using the asymptotic
    ! expressions by [Tow2003a]

    if (REAL(m*nu) >= 0._WP) then

       $ASSERT(k >= 0,Invalid k)

       if (k > 0) then

          s = k - 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))


       else

          lambda = m**2*(2._WP*m*nu)/(2._WP*m*nu - 1._WP)

       endif

    else

       if (k >= -1) then

          s = k + 1

          lambda = m*nu + m**2 + 0.5_WP*nu**2*(2*s + 1)**2* &
                                 (1._WP + SQRT(1._WP + 4._WP*(m*nu + m**2)/(nu**2*(2*s + 1)**2)))

!          lambda = m*nu + m**2 + 0.5_WP*(2*s + 1)**2* &
!                                 (nu**2 + SQRT(nu**4 + 4._WP*nu**2*(m*nu + m**2)/(2*s + 1)**2))

       else

          s = -k -1

          lambda = (m*nu - m**2)**2/(nu**2*(2*s+1)**2)

       endif

    endif

    ! Finish

    return

  end function lambda_asymp_c_

end module gyre_trad_func
