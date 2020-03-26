.. _osc-params:

Oscillation Parameters
======================

The :nml_g:`osc` namelist group defines oscillation parameters; the
input file can contain one or more, but only the last (matching) one
is used.  Allowable parameters are:

:nml_n:`inner_bound` (default :nml_v:`'REGULAR'`)
  Inner boundary conditions; one of:

  - :nml_v:`'REGULAR'` : Regularity-enforcing (only valid when inner grid point is at :math:`x = 0`)
  - :nml_v:`'ZERO_R'` : Zero radial displacement (only valid when inner grid point is at :math:`x \neq 0`)
  - :nml_v:`'ZERO_H'` : Zero horizontal displacement (only valid when inner grid point is at :math:`x \neq 0`)

:nml_n:`outer_bound` (default :nml_v:`'VACUUM'`)
  Outer boundary conditions; one of:

  - :nml_v:`'VACUUM'` : Zero surface pressure
  - :nml_v:`'DZIEM'` : Formulation following :ads_citet:`dziembowski:1971`
  - :nml_v:`'UNNO'` : Formulation following :ads_citet:`unno:1989`
  - :nml_v:`'JCD'` : Formulation following Jörgen Christensen-Dalsgaard (ADIPLS)

:nml_n:`variables_set` (default :nml_v:`'GYRE'`)
  Dependent variables in oscillation equations; one of:

  - :nml_v:`'GYRE'` : GYRE formulation, as desciribed in the :repo:`equations.pdf <doc/equations.pdf>` document
  - :nml_v:`'DZIEM'` : Formulation following :ads_citet:`dziembowski:1971`
  - :nml_v:`'JCD'` : Formulation following Jörgen Christensen-Dalsgaard (ADIPLS)
  - :nml_v:`'MIX'` : Mixed formulation (:nml_v:`'JCD'` for gravitational components, :nml_v:`'DZIEM'` for mechanical components)
  - :nml_v:`'LAGP'` : Lagrangian pressure perturbation formulation

:nml_n:`inertia_norm` (default :nml_v:`'BOTH'`)
  Inertia normalization factor; one of

  - :nml_v:`'RADIAL'` : Radial amplitude squared, :math:`|\xi_{\rm r}|^{2}`, evaluated at :nml_v:`x_ref`
  - :nml_v:`'HORIZ'` : Horizontal amplitude squared, :math:`|\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_v:`x_ref`
  - :nml_v:`'BOTH'` : Overall amplitude squared, :math:`|\xi_{\rm r}|^{2} + |\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_v:`x_ref`

:nml_n:`rotation_method` (default :nml_v:`'DOPPLER'`)
  rotation method; one of:

  - :nml_v:`'DOPPLER'` : Doppler shift
  - :nml_v:`'TAR'` : Traditional approximation of rotation

:nml_n:`time_factor` (default :nml_v:`'OSC'`)
  time-dependence factor in pulsation equations; one of:

  - :nml_v:`OSC` : Oscillatory, :math:`\propto \exp(-{\rm i} \omega t)`
  - :nml_v:`EXP` : Exponential, :math:`\propto \exp(-\omega t)`

:nml_n:`conv_scheme` (default :nml_v:`'FROZEN_PESNELL_1'``)
  convection treatment scheme; one of:

  - :nml_v:`'FROZEN_PESNELL_1'` : Freeze convective heating altogether;
    case 1 described by :ads_citet:`pesnell:1990`
  - :nml_v:`'FROZEN_PESNELL_4'` : Freeze Lagrangian perturbation of convective luminosity;
    case 4 described by :ads_citet:`pesnell:1990`

:nml_n:`deps_scheme` (default :nml_v:`'MODEL'`)
  scheme for calculating burning partial derivatives
  :math:`(\partial\ln\epsilon/\partial\ln T)_{\rho}` and
  :math:`(\partial\ln\epsilon/\partial\ln\rho)_{T}`; one of

  - :nml_v:`'MODEL'` : Use values from model
  - :nml_v:`'FILE'` : Use complex (phase-lagged) values from separate file

:nml_n:`deps_file` (default :nml_v:`''`)
  Name of epsilon partial derivatives file, when :nml_n:`deps_scheme` is :nml_v:`'FILE'`

:nml_n:`deps_file_format` (default :nml_v:`'WOLF'`)
  Format of epsilon partial derivative file, when :nml_n:`deps_scheme`
  is :nml_v:`'FILE'`; one of:

  - :nml_v:`'WOLF'` : Format used in preparation of :ads_citet:`wolf:2018`

:nml_n:`x_ref` (default :nml_v:`1` or outer grid point, whichever is smaller)
  Reference fractional radius for photosphere, normalizations etc.
   
:nml_n:`nonadiabatic` (default :nml_v:`.FALSE.`)
  Flag to include non-adiabatic effects
  
:nml_n:`quasiad_eigfuncs` (default :nml_v:`.FALSE.`)
  Flag to calculate quasi-adiabatic entropy/luminosity eigenfunctions
  during adiabatic calculations

:nml_n:`alpha_th` (defaualt :nml_v:`1.`)
  Scaling factor for thermal timescale in energy equation. Set to
  :nml_v:`0.` to recover the non-adiabatic reversible (NAR) limit, and
  to large values to approach the adiabatic limit

:nml_n:`cowling_approx` (default :nml_v:`.FALSE.`)
  Flag to use the Cowling approximation

:nml_n:`narf_approx` (default :nml_v:`.FALSE.`)
  Flag to use the non-adiabatic, radial flux (NARF) approximation
  
:nml_n:`eddington_approx` (default :nml_v:`.FALSE.`)
  Flag to use the Eddington approximation

:nml_n:`complex_lambda` (default :nml_v:`.FALSE.`)
  Flag to use complex arithmetic when evaluating angular eigenvalues
  lambda

:nml_n:`reduce_order` (default :nml_v:`.TRUE.`)
   Flag to reduce the order of the *adiabatic* radial-pulsation
   equations from 4 to 2

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
