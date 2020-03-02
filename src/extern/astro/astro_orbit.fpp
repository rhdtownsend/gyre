! Module   : astro_orbit
! Purpose  : binary orbits

$include 'core.inc'

module astro_orbit

  ! Uses

  use core_kinds
  use core_constants
  use core_func

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type orbit_t
     real(WP) :: M_pri
     real(WP) :: M_sec
     real(WP) :: period
     real(WP) :: e
   contains
     private
     procedure, public :: phase => phase_
     procedure, public :: a => a_
     procedure, public :: rho => rho_
     procedure, public :: rho_dot => rho_dot_
     procedure, public :: r => r_
     procedure, public :: r_dot => r_dot_
     procedure, public :: upsilon => upsilon_
     procedure, public :: upsilon_dot => upsilon_dot_
     procedure, public :: vr => vr_
     procedure, public :: vv => vv_
     procedure, public :: vr_L => vr_L_
  end type orbit_t

  type, extends (func_t) :: lagrange_func_t
     real(WP) :: M_pri
     real(WP) :: M_sec
     real(WP) :: Omega
     real(WP) :: x_pri
     real(WP) :: x_sec
   contains
     procedure, public :: eval_c_
  end type lagrange_func_t

  ! Acces specifiers

  private

  public :: orbit_t

  ! Procedures

contains

  function phase_ (this, t, t_peri) result (phase)

    class(orbit_t), intent(in) :: this
    real(WP), intent(in)       :: t
    real(WP), intent(in)       :: t_peri
    real(WP)                   :: phase

    ! Evaluate the orbital phase, relative to the given
    ! periastron time

    phase = TWOPI*MODULO((t - t_peri)/this%period, 1._WP)

    ! Finish

    return

  end function phase_

!****

  function a_ (this, star) result (a)

    class(orbit_t), intent(in)         :: this
    character(*), optional, intent(in) :: star
    real(WP)                           :: a

    ! Evaluate the semi-major axis

    a = ((this%period/TWOPI)**2*G_GRAVITY*(this%M_pri + this%M_sec))**(1._WP/3._WP)

    if (PRESENT(star)) then

       select case (star)
       case ('PRI')
          a = a*(this%M_sec/(this%M_pri + this%M_sec))
       case ('SEC')
          a = a*(this%M_pri/(this%M_pri + this%M_sec))
       case default
          $ABORT(Invalid star)
       end select

    end if

    ! Finish

    return

  end function a_

!****

  function rho_ (this, upsilon, star) result (rho)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: upsilon
    character(*), optional, intent(in) :: star
    real(WP)                           :: rho

    ! Evaluate the dimensionless radius

    if (PRESENT(star)) then

       select case (star)
       case ('PRI')
          rho = (1._WP - this%e**2)/(1._WP + this%e*COS(upsilon))
       case ('SEC')
          rho = (1._WP - this%e**2)/(1._WP - this%e*COS(upsilon))
       case default
          $ABORT(Invalid star)
       end select

    else

       rho = (1._WP - this%e**2)/(1._WP + this%e*COS(upsilon))

    endif

    ! Finish

    return

  end function rho_

!****

  function rho_dot_ (this, upsilon, star) result (rho_dot)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: upsilon
    character(*), optional, intent(in) :: star
    real(WP)                           :: rho_dot

    real(WP) :: rho

    ! Evaluate the rate of change of dimensionless 
    ! radius

    rho = this%rho(upsilon, star)

    rho_dot = (TWOPI/this%period)*SQRT(MAX(2._WP/rho - (1._WP - this%e**2)/rho**2 - 1._WP, 0._WP))

    ! Correct the sign

    if (PRESENT(star)) then

       select case (star)
       case ('PRI')
          if (upsilon < 0.) rho_dot = -rho_dot
       case ('SEC')
          if (upsilon > 0.) rho_dot = -rho_dot
       case default
          $ABORT(Invalid star)
       end select

    else

       if (upsilon < 0.) rho_dot = -rho_dot

    endif

    ! Finish

    return

  end function rho_dot_

!****

  function r_ (this, phase, star) result (r)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: r

    real(WP) :: upsilon

    ! Evaluate the radius

    upsilon = this%upsilon(phase, star)

    r = this%a(star)*this%rho(upsilon, star)

    ! Finish

    return

  end function r_

!****

  function r_dot_ (this, phase, star) result (r_dot)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: r_dot

    real(WP) :: upsilon

    ! Evaluate the rate of change of radius

    upsilon = this%upsilon(phase, star)

    r_dot = this%a(star)*this%rho_dot(upsilon, star)

    ! Finish

    return

  end function r_dot_

!****

  function upsilon_ (this, phase, star) result (upsilon)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: upsilon

    real(WP) :: M
    real(WP) :: E

    ! Evaluate the true anomaly

    M = ATAN2(SIN(phase), COS(phase))

    E = E_5_(REAL(this%e, DP), REAL(M, DP))

    upsilon = 2._WP*ATAN(SQRT((1._WP + this%e)/(1._WP - this%e))*TAN(E/2._WP))

    if (PRESENT(star)) then

       select case (star)
       case ('PRI')
       case ('SEC')
          upsilon = ATAN2(-SIN(upsilon), -COS(upsilon))
       case default
          $ABORT(Invalid star)
       end select

    endif

    ! Finish

    return

  contains

    function E_5_ (e, M) result (E_5)

      real(DP), intent(in) :: e
      real(DP), intent(in) :: M
      real(DP)             :: E_5

      real(DP) :: alpha
      real(DP) :: d
      real(DP) :: q
      real(DP) :: r
      real(DP) :: w
      real(DP) :: E_1
      real(DP) :: M_ast
      real(DP) :: f_0
      real(DP) :: f_1
      real(DP) :: f_2
      real(DP) :: f_3
      real(DP) :: f_4
      real(DP) :: delta_3
      real(DP) :: delta_4
      real(DP) :: delta_5

      ! Given the eccentriciy and the mean anomaly, find the eccentric
      ! anomaly using the fifth- order approximation by Landis Markley
      ! (1995, Celestial Mechanics and Dynamical Astronomy, 63,
      ! 101-111).

      ! Calculate the Pade alpha parameter

      alpha = (3._DP*PI**2 + 1.6_DP*PI*(PI - ABS(M))/(1._DP + e))/ &
              (PI**2 - 6._DP)

      ! Set up the first approximation to the eccentric 
      ! anomaly (LM, eqn. 15)

      d = 3._DP*(1._DP - e) + alpha*e

      q = 2._DP*alpha*d*(1._DP - e) - M**2
      r = 3._DP*alpha*d*(d - 1._DP + e)*M + M**3

      w = (ABS(r) + SQRT(q**3 + r**2))**(2._DP/3._DP)

      E_1 = (2._DP*r*w/(w**2 + w*q + q**2) + M)/d

      ! Refine this approximation using the fifth-order correction

      if (this%e > 0.5_DP .AND. E_1 < 1._DP) then
         M_ast = (1._DP - this%e)*E_1 + this%e*E_1**3*(      &
                 -1.7454287843856404E-6_DP*E_1**6 +          &
                  4.1584640418181644E-4_DP*E_1**4 -          &
                  3.0956446448551138E-2_DP*E_1**2 + 1._DP)/( &
                  1.7804367119519884E-8_DP*E_1**8 +          &
                  5.9727613731070647E-6_DP*E_1**6 +          &
                  1.0652873476684142E-3_DP*E_1**4 +          &
                  1.1426132130869317E-1_DP*E_1**2 + 6._DP)
      else
         M_ast = E_1 - this%e*SIN(E_1)
      endif

      f_0 = M_ast - M
      f_1 = 1._DP - this%e*COS(E_1)
      f_2 = E_1 - M_ast
      f_3 = 1._DP - f_1
      f_4 = -f_2

      delta_3 = -f_0/(f_1 - f_0*f_2/(2._DP*f_1))
      delta_4 = -f_0/(f_1 + delta_3*f_2/2._DP + delta_3**2*f_3/6._DP)
      delta_5 = -f_0/(f_1 + delta_4*f_2/2._DP + delta_4**2*f_3/6._DP + delta_4**3*f_4/24._DP)

      E_5 = E_1 + delta_5

      ! Finish

      return

    end function E_5_

  end function upsilon_

!****

  function upsilon_dot_ (this, phase, star) result (upsilon_dot)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: upsilon_dot

    real(WP) :: upsilon

    ! Evaluate the rate-of-change of the true anomaly

    upsilon = this%upsilon(phase, star)

    upsilon_dot = (TWOPI/this%period)*SQRT(1._WP - this%e**2)/this%rho(upsilon, star)**2

    ! Finish

    return

  end function upsilon_dot_

!****

  function vr_ (this, phase, star) result (vr)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: vr(3)

    real(WP) :: upsilon
    real(WP) :: r

    ! Evaluate the Cartesian position vector

    upsilon = this%upsilon(phase, star)
    r = this%a(star)*this%rho(upsilon, star)

    vr = r*[COS(upsilon),SIN(upsilon),0._WP]

    ! Finish

    return

  end function vr_

!****

  function vv_ (this, phase, star) result (vv)

    class(orbit_t), intent(in)         :: this
    real(WP), intent(in)               :: phase
    character(*), optional, intent(in) :: star
    real(WP)                           :: vv(3)

    real(WP) :: upsilon
    real(WP) :: r
    real(WP) :: v_r
    real(WP) :: v_t

    ! Evaluate the Cartesian velocity vector

    upsilon = this%upsilon(phase, star)
    r = this%a(star)*this%rho(upsilon, star)

    v_r = this%r_dot(phase, star)
    v_t = r*this%upsilon_dot(phase, star)

    vv = [v_r*COS(upsilon) - v_t*SIN(upsilon), &
          v_r*SIN(upsilon) + v_t*COS(upsilon), &
          0._WP]

    ! Finish

    return

  end function vv_

!****

  function vr_L_ (this, phase, i) result (vr_L)

    class(orbit_t), intent(in) :: this
    real(WP), intent(in)       :: phase
    integer, intent(in)        :: i
    real(WP)                   :: vr_L(3)

    real(WP)              :: upsilon
    real(WP)              :: rho
    real(WP)              :: r_pri
    real(WP)              :: r_sec
    type(lagrange_func_t) :: lf
    real(WP)              :: x_a
    real(WP)              :: x_b
    real(WP)              :: f_a
    real(WP)              :: f_b
    real(WP)              :: x
    real(WP)              :: r
    real(WP)              :: delta_upsilon

    ! Calculate the Cartesian position vector of the i'th
    ! Lagrangian point

    upsilon = this%upsilon(phase)
    rho = this%rho(upsilon)

    r_pri = this%a(star='PRI')*rho
    r_sec = this%a(star='SEC')*rho

    select case(i)

    case(1,2,3)

       ! In-line points

       ! Set up the lagrange function module variable

       lf%M_pri = this%M_pri
       lf%M_sec = this%M_sec

       lf%Omega = this%upsilon_dot(phase)

       lf%x_pri = r_pri
       lf%x_sec = -r_sec

       ! Locate the point by finding the appropriate root of the force
       ! function

       select case(i)
       case(1)

          x_a = lf%x_pri
          x_b = lf%x_sec

          f_a = HUGE(0._WP)
          f_b = -HUGE(0._WP)

       case(2)

          x_a = lf%x_sec
          x_b = lf%x_sec - (G_GRAVITY*(lf%M_pri+lf%M_sec))**(1._WP/3._WP)*lf%Omega**(-2._WP/3._WP)

          f_a = HUGE(0._WP)
          f_b = lf%eval(x_b)

       case(3)

          x_a = lf%x_pri
          x_b = lf%x_pri + (G_GRAVITY*(lf%M_pri+lf%M_sec))**(1._WP/3._WP)*lf%Omega**(-2._WP/3._WP)

          f_a = -HUGE(0._WP)
          f_b = lf%eval(x_b)

       end select

       x = lf%root(x_a, x_b, 0._WP, f_a, f_b)

       ! Calculate the Cartesian coordinates of the point

       vr_L = x*[COS(upsilon),SIN(upsilon),0._WP]

    case(4,5)

       ! Trojan points

       ! Set up the radius and offset angle of the vertex of the
       ! equilateral triangle

       r = SQRT(r_pri**2 + (r_pri + r_sec)**2 - r_pri*(r_pri + r_sec))
       delta_upsilon = ACOS((r**2 + r_pri**2 - (r_pri + r_sec)**2)/(2._WP*r*r_pri))

       ! Calculate the Cartesian coordinates of the point

       select case(i)
       case(4)
          vr_L = r*[COS(upsilon + delta_upsilon),SIN(upsilon + delta_upsilon),0._WP]
       case(5)
          vr_L = r*[COS(upsilon - delta_upsilon),SIN(upsilon - delta_upsilon),0._WP]
       end select
       
    case default

       $ABORT_DEBUG(Invalid i)

    end select

    ! Finish

    return

  end function vr_L_

!****

  function eval_c_ (this, z) result (f_z)

    class(lagrange_func_t), intent(inout) :: this
    complex(WP), intent(in)               :: z
    complex(WP)                           :: f_z

    ! Evaluate the function for lagrange_func_t

    f_z = lagrange_force(this, REAL(z))

    ! Finish

    return

  end function eval_c_

!****

  function lagrange_force (lf, x) result (f)

    class(lagrange_func_t), intent(in) :: lf
    real(WP), intent(in)               :: x
    real(WP)                           :: f

    ! Calculate the force on a particle at position x

    if (x < lf%x_sec) then

       f = G_GRAVITY*lf%M_pri/(x - lf%x_pri)**2 + G_GRAVITY*lf%M_sec/(x - lf%x_sec)**2

    elseif (x > lf%x_pri) then

       f = -G_GRAVITY*lf%M_pri/(x - lf%x_pri)**2 - G_GRAVITY*lf%M_sec/(x - lf%x_sec)**2

    else

       f = G_GRAVITY*lf%M_pri/(x - lf%x_pri)**2 - G_GRAVITY*lf%M_sec/(x - lf%x_sec)**2

    end if

    f  = f + lf%Omega**2*x

    ! Finish

    return

  end function lagrange_force

end module astro_orbit
