! Program  : gyre_trad_fit
! Purpose  : fits to traditional approximation eigenvalues
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

module gyre_trad_fit

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel

  use gyre_trad_eigen
  use gyre_cheb_fit

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: trad_fit_t
     private
     type(cheb_fit_t) :: cf_neg
     type(cheb_fit_t) :: cf_ctr
     type(cheb_fit_t) :: cf_pos
     logical          :: df_neg
     logical          :: df_ctr
     logical          :: df_pos
     real(WP)         :: nu_t
     integer, public  :: m
     integer, public  :: k
   contains
     private
     procedure       :: lambda_r_
     procedure       :: lambda_c_
     generic, public :: lambda => lambda_r_, lambda_c_
  end type trad_fit_t

  ! Interfaces

  interface trad_fit_t
     module procedure trad_fit_t_tol_
  end interface trad_fit_t

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

  public :: trad_fit_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif

  ! Procedures

contains

  function trad_fit_t_tol_ (m, k, cheb_tol) result (tf)

    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: cheb_tol
    type(trad_fit_t)     :: tf

    integer :: l

    ! Construct the trad_fit_t with the specified tolerances

    tf%m = m
    tf%k = k

    if (k >= 0) then

       ! Gravito-inertial waves

       tf%nu_t = 1._WP

       tf%cf_neg = cheb_fit_t(-1._WP, 0._WP, cheb_tol, f_grav_outer_)
       tf%cf_ctr = cheb_fit_t(-1._WP, 1._WP, cheb_tol, f_grav_inner_)
       tf%cf_pos = cheb_fit_t( 0._WP, 1._WP, cheb_tol, f_grav_outer_)

       tf%df_neg = .TRUE.
       tf%df_ctr = .TRUE.
       tf%df_pos = .TRUE.

    else

       ! Rossby waves

       l = ABS(m) + ABS(k) - 1

       if (m > 0) then

          tf%nu_t = REAL(l*(l+1), WP)/REAL(m, WP)

          tf%cf_neg = cheb_fit_t(-1._WP, 0._WP, cheb_tol, f_ross_)

          tf%df_neg = .TRUE.
          tf%df_pos = .FALSE.

       elseif (m < 0) then

          tf%nu_t = REAL(l*(l+1), WP)/REAL(-m, WP)

          tf%cf_pos = cheb_fit_t(0._WP, 1._WP, cheb_tol, f_ross_)

          tf%df_neg = .FALSE.
          tf%df_pos = .TRUE.

       else

          tf%df_neg = .FALSE.
          tf%df_pos = .FALSE.

       endif

       tf%df_ctr = .FALSE.

    endif
       
    ! Finish

    return

  contains

    function f_grav_inner_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the gravity-wave eigenvalue function in the inner
      ! region (|nu| < nu_t)

      if (m == 0._WP .AND. k == 0._WP) then
         
         f = 1._WP

      else

         l = ABS(m) + ABS(k) - 1

         associate (nu => tf%nu_t*x)
           f = lambda(nu, m, k)/lambda_norm_grav_inner_(nu, m, k)
         end associate

      endif

      ! Finish

      return

    end function f_grav_inner_

    function f_grav_outer_ (x) result (f)

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
              f = lambda(nu, m, k)/lambda_norm_grav_outer_(nu, m, k)
            end associate

         endif

      endif

      ! Finish

      return

    end function f_grav_outer_

    function f_ross_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the Rossby-wave eigenvaule function

      if (x == 0._WP) then

         f = 1._WP

      elseif (m > 0 .AND. x == -1._WP) then

         f = 0._WP

      elseif (m < 0 .AND. x == 1._WP) then

         f = 0._WP

      else

         associate (nu => tf%nu_t/x)

           if (k == -1) then
              f = lambda(nu, m, k)/lambda_norm_grav_outer_(nu, m, k)
           else
              f = lambda(nu, m, k)/lambda_norm_ross_(nu, m, k)
           endif

         end associate

      endif

    end function f_ross_

  end function trad_fit_t_tol_

  !****

  $if ($HDF5)

  subroutine read_ (hg, tf)

    type(hgroup_t), intent(inout) :: hg
    type(trad_fit_t), intent(out) :: tf

    type(hgroup_t) :: hg_comp

    ! Read the trad_fit_t

    call read_attr(hg, 'm', tf%m)
    call read_attr(hg, 'k', tf%k)

    call read_attr(hg, 'nu_t', tf%nu_t)

    call read_attr(hg, 'df_neg', tf%df_neg)
    call read_attr(hg, 'df_ctr', tf%df_ctr)
    call read_attr(hg, 'df_pos', tf%df_pos)

    if (tf%df_neg) then
       hg_comp = hgroup_t(hg, 'cf_neg')
       call read(hg_comp, tf%cf_neg)
       call hg_comp%final()
    endif

    if (tf%df_ctr) then
       hg_comp = hgroup_t(hg, 'cf_ctr')
       call read(hg_comp, tf%cf_ctr)
       call hg_comp%final()
    endif

    if (tf%df_pos) then
       hg_comp = hgroup_t(hg, 'cf_pos')
       call read(hg_comp, tf%cf_pos)
       call hg_comp%final()
    endif

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, tf)

    type(hgroup_t), intent(inout) :: hg
    type(trad_fit_t), intent(in)  :: tf

    type(hgroup_t) :: hg_comp

    ! Write the trad_fit_t

    call write_attr(hg, 'm', tf%m)
    call write_attr(hg, 'k', tf%k)

    call write_attr(hg, 'nu_t', tf%nu_t)

    call write_attr(hg, 'df_neg', tf%df_neg)
    call write_attr(hg, 'df_ctr', tf%df_ctr)
    call write_attr(hg, 'df_pos', tf%df_pos)

    if (tf%df_neg) then
       hg_comp = hgroup_t(hg, 'cf_neg')
       call write(hg_comp, tf%cf_neg)
       call hg_comp%final()
    endif

    if (tf%df_ctr) then
       hg_comp = hgroup_t(hg, 'cf_ctr')
       call write(hg_comp, tf%cf_ctr)
       call hg_comp%final()
    endif

    if (tf%df_pos) then
       hg_comp = hgroup_t(hg, 'cf_pos')
       call write(hg_comp, tf%cf_pos)
       call hg_comp%final()
    endif

    ! Finish

    return

  end subroutine write_

  $endif

  !****

  function lambda_r_ (this, nu) result (lambda)

    class(trad_fit_t), intent(in), target :: this
    real(WP), intent(in)                  :: nu
    real(WP)                              :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (real)

    if (this%k >= 0) then

       ! Gravity waves

       if (nu <= -this%nu_t) then
          lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_grav_outer_(nu, this%m, this%k)
       elseif (nu >= this%nu_t) then
          lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_grav_outer_(nu, this%m, this%k)
       else
          lambda = this%cf_ctr%eval(nu/this%nu_t)*lambda_norm_grav_inner_(nu, this%m, this%k)
       endif

    else

       ! Rossby waves

       if (this%m > 0) then

          $ASSERT(nu <= -this%nu_t,Invalid nu for Rossby waves)

          if (this%k == -1) then
             lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_grav_outer_(nu, this%m, this%k)
          else
             lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_ross_(nu, this%m, this%k)
          endif

       elseif (this%m < 0) then

          $ASSERT(nu >= this%nu_t,Invalid nu for Rossby waves)

          if (this%k == -1) then
             lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_grav_outer_(nu, this%m, this%k)
          else
             lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_ross_(nu, this%m, this%k)
          endif

       else

          $ABORT(Invalid m for Rossby waves)

       endif

    end if

    ! Finish

    return

  end function lambda_r_

  !****

  function lambda_c_ (this, nu) result (lambda)

    class(trad_fit_t), intent(in) :: this
    complex(WP), intent(in)        :: nu
    complex(WP)                    :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (complex)

    if (this%k >= 0) then

       ! Gravity waves

       if (REAL(nu) < -this%nu_t) then
          lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_grav_outer_(REAL(nu), this%m, this%k)
       elseif (REAL(nu) > this%nu_t) then
          lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_grav_outer_(REAL(nu), this%m, this%k)
       else
          lambda = this%cf_ctr%eval(nu/this%nu_t)*lambda_norm_grav_inner_(REAL(nu), this%m, this%k)
       endif

    else

       ! Rossby waves

       if (this%m > 0) then

          $ASSERT(REAL(nu) <= -this%nu_t,Invalid nu for Rossby waves)

          if (this%k == -1) then
             lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_grav_outer_(REAL(nu), this%m, this%k)
          else
             lambda = this%cf_neg%eval(this%nu_t/nu)*lambda_norm_ross_(REAL(nu), this%m, this%k)
          endif

       else

          $ASSERT(REAL(nu) >= this%nu_t,Invalid nu for Rossby waves)

          if (this%k == -1) then
             lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_grav_outer_(REAL(nu), this%m, this%k)
          else
             lambda = this%cf_pos%eval(this%nu_t/nu)*lambda_norm_ross_(REAL(nu), this%m, this%k)
          endif

       end if

    end if

    ! Finish

    return

  end function lambda_c_

  !****

  function lambda_norm_grav_inner_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    integer :: l

    ! Evaluate the gravity-wave eigenvalue normalization function in
    ! the inner region (|nu| < nu_t)

    l = ABS(m) + k

    lambda_norm = l*(l+1)

    ! Finish

    return
    
  end function lambda_norm_grav_inner_
  
  !****

  function lambda_norm_grav_outer_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    integer :: s

    ! Evaluate the gravity-wave eigenvalue normalization function in
    ! the outer region (|nu| > nu_t)

    if (m*nu >= 0._WP) then

       if (k > 0) then
          s = k - 1
          lambda_norm = nu**2*(2*s + 1)**2
       else
          lambda_norm = m**2
       endif

    else

       s = k + 1
       lambda_norm = nu**2*(2*s + 1)**2

    endif

    ! Finish

    return
    
  end function lambda_norm_grav_outer_
  
  !****

  function lambda_norm_ross_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    integer :: s

    ! Evaluate the Rossby-wave eigenvalue normalization function

    s = -k -1

    lambda_norm = REAL(m, WP)**2/(2*s+1)**2

    ! Finish

    return
    
  end function lambda_norm_ross_

end module gyre_trad_fit
