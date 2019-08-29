Oscillation Parameters
======================

The :nml_g:`osc` namelist group defines oscillation parameters; the
input file can contain one or more, but only the last (matching) one
is used.  Allowable fields are:

:nml_o:`inner_bound` (default :nml_l:`'REGULAR'`)
  Inner boundary conditions; one of:

  - :nml_l:`'REGULAR'` : Regularity-enforcing (only valid when inner grid point is at :math:`x = 0`)
  - :nml_l:`'ZERO_R'` : Zero radial displacement (only valid when inner grid point is at :math:`x \neq 0`)
  - :nml_l:`'ZERO_H'` : Zero horizontal displacement (only valid when inner grid point is at :math:`x \neq 0`)

:nml_o:`outer_bound` (default :nml_l:`'VACUUM'`)
  Outer boundary conditions; one of:

  - :nml_l:`'VACUUM'` : Zero surface pressure
  - :nml_l:`'DZIEM'` : Formulation following :cite:`Dziembowski:1971a`
  - :nml_l:`'UNNO'` : Formulation following :cite:`Unno:1989a`
  - :nml_l:`'JCD'` : Formulation following Jörgen Christensen-Dalsgaard (ADIPLS)

:nml_o:`variables_set` (default :nml_l:`'GYRE'`)
  Dependent variables in oscillation equations; one of:

  - :nml_l:`'GYRE'` : GYRE formulation, as desciribed in the :repo:`equations.pdf <doc/equations.pdf>` document
  - :nml_l:`'DZIEM'` : Formulation following :cite:`Dziembowski:1971a`
  - :nml_l:`'JCD'` : Formulation following Jörgen Christensen-Dalsgaard (ADIPLS)
  - :nml_l:`'MIX'` : Mixed formulation (:nml_l:`'JCD'` for gravitational components, :nml_l:`'DZIEM'` for mechanical components)
  - :nml_l:`'LAGP'` : Lagrangian pressure perturbation formulation

:nml_o:`inertia_norm` (default :nml_l:`'BOTH'`)
  Inertia normalization factor; one of

  - :nml_l:`'RADIAL'` : Radial amplitude squared, :math:`|\xi_{\rm r}|^{2}`, evaluated at :nml_l:`x_ref`
  - :nml_l:`'HORIZ'` : Horizontal amplitude squared, :math:`|\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_l:`x_ref`
  - :nml_l:`'BOTH'` : Overall amplitude squared, :math:`|\xi_{\rm r}|^{2} + |\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_l:`x_ref`

:nml_o:`rotation_method` (default :nml_l:`'DOPPLER'`)
  rotation method; one of:

  - :nml_l:`'DOPPLER'` : Doppler shift
  - :nml_l:`'TAR'` : Traditional approximation of rotation

:nml_o:`time_factor` (default :nml_l:`'OSC'`)
  time-dependence factor in pulsation equations; one of:

  - :nml_l:`OSC` : Oscillatory, :math:`\propto \exp(-{\rm i} \omega t)`
  - :nml_l:`EXP` : Exponential, :math:`\propto \exp(-\omega t)`

:nml_o:`conv_scheme` (default :nml_l:`'FROZEN_PESNELL_1'``)
  convection treatment scheme; one of:

  - :nml_l:`'FROZEN_PESNELL_1'` : Freeze convective heating altogether;
    case 1 described by :cite:`Pesnell:1990a`
  - :nml_l:`'FROZEN_PESNELL_4'` : Freeze Lagrangian perturbation of convective luminosity;
    case 4 described by :cite:`Pesnell:1990a`

:nml_o:`deps_scheme` (default :nml_l:`'MODEL'`)
  scheme for calculating burning partial derivatives
  :math:`(\partial\ln\epsilon/\partial\ln T)_{\rho}` and
  :math:`(\partial\ln\epsilon/\partial\ln\rho)_{T}`; one of

  - :nml_l:`'MODEL'` : Use values from model
  - :nml_l:`'FILE'` : Use complex (phase-lagged) values from separate file

:nml_o:`deps_file` (default :nml_l:`''`)
  Name of epsilon partial derivatives file, when :nml_o:`deps_scheme` is :nml_l:`'FILE'`

:nml_o:`deps_file_format` (default :nml_l:`'WOLF'`)
  Format of epsilon partial derivative file, when :nml_o:`deps_scheme`
  is :nml_l:`'FILE'`; one of:

  - :nml_l:`'WOLF'` : Format used in preparation of :cite:`Wolf:2018a`

:nml_o:`x_ref` (default :nml_l:`1` or outer grid point, whichever is smaller)
  Reference fractional radius for photosphere, normalizations etc.
   
:nml_o:`nonadiabatic` (default :nml_l:`.FALSE.`)
  Flag to include non-adiabatic effects
  
:nml_o:`quasiad_eigfuncs` (default :nml_l:`.FALSE.`)
  Flag to calculate quasi-adiabatic entropy/luminosity eigenfunctions
  during adiabatic calculations

:nml_o:`cowling_approx` (default :nml_l:`.FALSE.`)
  Flag to use the Cowling approximation

:nml_o:`nar_approx` (default :nml_l:`.FALSE.`)
  Flag to use the non-adiabatic reversible (NAR) approximation
  
:nml_o:`narf_approx` (default :nml_l:`.FALSE.`)
  Flag to use the non-adiabatic, radial flux (NARF) approximation
  
:nml_o:`eddington_approx` (default :nml_l:`.FALSE.`)
  Flag to use the Eddington approximation

:nml_o:`complex_lambda` (default :nml_l:`.FALSE.`)
  Flag to use complex arithmetic when evaluating angular eigenvalues
  lambda

:nml_o:`reduce_order` (default :nml_l:`.TRUE.`)
   Flag to reduce the order of the *adiabatic* radial-pulsation
   equations from 4 to 2

:nml_o:`tag_list` (default :nml_l:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
