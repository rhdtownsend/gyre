! Program : eval_hough
! Purpose : evaluate Hough functions

$include 'core.inc'

program eval_hough

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup
  use core_system

  use astro_hough

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  real, parameter :: LAMBDA_TOL = 1E-5

  ! Variables

  character(:), allocatable :: filename
  integer                   :: m
  integer                   :: k
  integer                   :: n_mu
  integer                   :: n_phi
  real(WP)                  :: nu
  integer                   :: i_mu
  integer                   :: i_phi
  real(WP), allocatable     :: mu(:)
  real(WP), allocatable     :: phi(:)
  integer                   :: i
  type(hough_t)             :: ho
  real(WP), allocatable     :: Theta(:,:)
  complex(WP), allocatable  :: xi(:,:,:)
  type(hgroup_t)            :: hg

  ! Read command-line arguments

  $ASSERT(n_arg() == 6,Syntax: eval_hough filename m k nu n_mu n_phi)

  call get_arg(1, filename)
  call get_arg(2, m)
  call get_arg(3, k)
  call get_arg(4, nu)
  call get_arg(5, n_mu)
  call get_arg(6, n_phi)

  ! Set up the mu and phi axes
  
  allocate(mu(n_mu))
  allocate(phi(n_phi))

  mu = [((-(n_mu-i) + (i-1))/(n_mu-1.),i=1,n_mu)]
  phi = [((-PI*(n_phi-i) + PI*(i-1))/(n_phi-1),i=1,n_phi)]

  ! Create the hough_t

  ho = hough_t(nu, m, k, 100)!LAMBDA_TOL)

  ! Calculate the displacement functions

  allocate(Theta(3,n_mu))
  allocate(xi(3,n_mu,n_phi))

  mu_loop : do i_mu = 1, n_mu
     Theta(:,i_mu) = ho%Theta(mu(i_mu))
     phi_loop : do i_phi = 1, n_phi
        xi(:,i_mu,i_phi) = ho%Theta_xi(mu(i_mu))*[CMPLX( 1._WP,  0._WP, WP), &
                                                  CMPLX( 1._WP,  0._WP, WP), &
                                                  CMPLX( 0._WP, -1._WP, WP)]*EXP(CMPLX(0._WP, 1._WP, KIND=WP)*m*phi(i_phi))
     end do phi_loop
  end do mu_loop

  ! Write out the results

  hg = hgroup_t(filename, CREATE_FILE)

  call write_attr(hg, 'm', m)
  call write_attr(hg, 'k', k)
  call write_attr(hg, 'nu', nu)

  call write_attr(hg, 'n_mu', n_mu)
  call write_attr(hg, 'n_phi', n_phi)

  call write_attr(hg, 'lambda', ho%lambda)

  call write_dset(hg, 'mu', mu)
  call write_dset(hg, 'phi', phi)

  call write_dset(hg, 'Theta_r', Theta(1,:))
  call write_dset(hg, 'Theta_t', Theta(2,:))
  call write_dset(hg, 'Theta_p', Theta(3,:))

  call write_dset(hg, 'xi_r', xi(1,:,:))
  call write_dset(hg, 'xi_t', xi(2,:,:))
  call write_dset(hg, 'xi_p', xi(3,:,:))

  call hg%final()

  ! Finish

end program eval_hough
