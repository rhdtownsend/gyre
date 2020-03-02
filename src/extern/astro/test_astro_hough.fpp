! Program : test_astro_hough
! Purpose : test harness for astro_hough module

$include 'core.inc'

program test_astro_hough

  ! Uses

  use core_kinds
  use core_spline

  use astro_hough

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: M_MIN = -5
  integer, parameter :: M_MAX = 5
  integer, parameter :: K_MIN = 0
  integer, parameter :: K_MAX = 5

  real(WP), parameter :: TOL = 1E-4_WP

  real(WP), parameter :: NU_MIN = -10._WP
  real(WP), parameter :: NU_MAX = 10._WP
  integer, parameter  :: N_NU = 100

  integer, parameter :: N_MU = 1000

  ! Variables

  integer       :: m
  integer       :: k
  integer       :: i_nu
  integer       :: i_mu
  real(WP)      :: nu
  integer       :: l_j
  integer       :: j
  integer       :: j_res
  real(WP)      :: mu(N_MU)
  type(hough_t) :: ho
  real(WP)      :: Theta(3)
  real(WP)      :: Theta_xi(3)
  real(WP)      :: Thetar(N_MU)
  real(WP)      :: Thetat(N_MU)
  real(WP)      :: Thetap(N_MU)
  real(WP)      :: Thetar_xi(N_MU)
  real(WP)      :: Thetat_xi(N_MU)
  real(WP)      :: Thetap_xi(N_MU)
  type(spline_t) :: sp_Thetar
  type(spline_t) :: sp_Thetat
  type(spline_t) :: sp_Thetap
  real(WP)      :: dThetar(N_MU)
  real(WP)      :: dThetat(N_MU)
  real(WP)      :: dThetap(N_MU)
  real(WP)      :: res_17(N_MU)
  real(WP)      :: res_18(N_MU)
  real(WP)      :: res_19(N_MU)

  ! First calculate the hough functions

  mu = [((-1._WP*(N_MU-i_mu) + 1._WP*(i_mu-1))/(N_MU-1),i_mu=1,N_MU)]

  ! Loop over m, k & nu

  do m = M_MIN, M_MAX
     do k = K_MIN, K_MAX
        do i_nu = 1, N_NU

           nu = (NU_MIN*(N_NU-i_nu) + NU_MAX*(i_nu-1))/(N_NU-1)

           ! Test for rossby-mode resonance

           if (MOD(k, 2) == 0) then
              l_j = ABS(m)
           else
              l_j = ABS(m) + 1
           endif

           j_res = 0

           scan_loop : do j = 1, 500
              if ((l_j+1)*(l_j+2) == -m*nu) then
                 j_res = j
                 exit scan_loop
              endif
              l_j = l_j + 2
           end do scan_loop

           if (j_res /= 0) then
              nu = nu*(1._WP + SQRT(EPSILON(0._WP)))
              print *,'Nudged'
           endif

           ! Construct the hough_t

           ho = hough_t(nu, m, k, tol)

           ! Loop over nu

           do i_mu = 1, N_MU

              ! Evaluate the Hough functions

              Theta = ho%Theta(mu(i_mu))
              Theta_xi = ho%Theta_xi(mu(i_mu))

              Thetar(i_mu) = Theta(1)
              Thetat(i_mu) = Theta(2)
              Thetap(i_mu) = Theta(3)

              Thetar_xi(i_mu) = Theta_xi(1)
              Thetat_xi(i_mu) = Theta_xi(2)
              Thetap_xi(i_mu) = Theta_xi(3)

           end do

           ! Check the extent to which the Hough functions satisfy Laplace's tidal equation

           ! First, generate spline fits

           sp_Thetar = spline_t(mu, Thetar, 'MONO')
           sp_Thetat = spline_t(mu, Thetat, 'MONO')
           sp_Thetap = spline_t(mu, Thetap, 'MONO')

           ! Calculate derivatives

           DThetar = (1._WP-mu**2)*sp_Thetar%deriv()
           DThetat = (1._WP-mu**2)*sp_Thetat%deriv()
           DThetap = (1._WP-mu**2)*sp_Thetap%deriv()

           ! Calculate the residuals of the tidal eqns. (17-19) in Townsend (2003)

           res_17 = -Thetat - nu*mu*Thetap - DThetar
           res_18 = -Thetap - nu*mu*Thetat - m*Thetar
           res_19 = ho%lambda*(1._WP-mu**2)*Thetar - DThetat + m*Thetap

           print *,'Max dev:', MAXVAL(ABS(res_17)), MAXVAL(ABS(res_18)), MAXVAL(ABS(res_19))

        end do
     end do
  end do

  ! Finish

end program test_astro_hough


              

  

  
