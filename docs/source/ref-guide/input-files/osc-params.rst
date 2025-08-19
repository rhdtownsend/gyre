.. _osc-params:

Oscillation Parameters
======================

The :nml_g:`osc` namelist group defines oscillation parameters, as follows:

:nml_n:`inner_bound` (default :nml_v:`'REGULAR'`)
  Inner boundary conditions; one of:

  - :nml_v:`'REGULAR'` : Regularity-enforcing (only valid when inner grid point is at :math:`x = 0`)
  - :nml_v:`'ZERO_R'` : Zero radial displacement (only valid when inner grid point is at :math:`x \ne 0`)
  - :nml_v:`'ZERO_H'` : Zero horizontal displacement (only valid when inner grid point is at :math:`x \ne 0`)

:nml_n:`outer_bound` (default :nml_v:`'VACUUM'`)
  Outer boundary conditions; one of:

  - :nml_v:`'VACUUM'` : Vanishing surface density
  - :nml_v:`'DZIEM'` : Formulation following :ads_citet:`dziembowski:1971`
  - :nml_v:`'UNNO'` : Formulation following :ads_citet:`unno:1989`
  - :nml_v:`'JCD'` : Formulation following Jørgen Christensen-Dalsgaard (ADIPLS)
  - :nml_v:`'ISOTHERMAL'` : Formulation based on local dispersion analysis for isothermal atmosphere
  - :nml_v:`'GAMMA'` : Vanishing displacement and derivative at outer boundary, intended for use with :math:`\gamma` modes (isolated g modes; see :ads_citealp:`ong:2020`)

:nml_n:`outer_bound_cutoff` (default :nml_v:`''`)
  Outer boundary conditions to use when evaluating cutoff frequencies (see :nml_n:`freq_units`); same options
  as :nml_n:`outer_bound`, and if left blank then takes its value from :nml_n:`outer_bound`

:nml_n:`outer_bound_branch` (default :nml_v:`'E_NEG'`)
  Dispersion relation solution branch to use for outer boundary
  conditions (when :nml_n:`outer_bound`\ =\ :nml_v:`'UNNO'`\ \|\ :nml_v:`'JCD'`\ \|\ :nml_v:`'ISOTHERMAL'`);
  one of

  - :nml_v:`'E_NEG'` : Outward-decaying energy density
  - :nml_v:`'E_POS'` : Outward-growing energy density
  - :nml_v:`'F_NEG'` : Outward energy flux
  - :nml_v:`'F_POS'` : Inward energy flux
  - :nml_v:`'V_NEG'` : Outward phase velocity
  - :nml_v:`'V_POS'` : Inward phase velocity

:nml_n:`variables_set` (default :nml_v:`'GYRE'`)
  Dependent variables in oscillation equations; one of:

  - :nml_v:`'GYRE'` : GYRE formulation, as described in the :ref:`osc-dimless-form` section
  - :nml_v:`'DZIEM'` : Formulation following :ads_citet:`dziembowski:1971`
  - :nml_v:`'JCD'` : Formulation following Jørgen Christensen-Dalsgaard (ADIPLS)
  - :nml_v:`'MIX'` : Mixed formulation (:nml_v:`'JCD'` for :math:`y_{3,4}`, :nml_v:`'DZIEM'` for :math:`y_{1,2}`)
  - :nml_v:`'LAGP'` : Lagrangian pressure perturbation formulation

:nml_n:`alpha_grv` (default :nml_v:`1.`)
  Scaling factor for gravitational potential perturbations (see the :math:`\alphagrv`
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_gbc` (default :nml_v:`1.`)
  Scaling factor for the displacement term in the outer gravitational potential boundary
  condition (see the :math:`\alphagbc` entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_thm` (default :nml_v:`1.`)
  Scaling factor for the thermal timescale (see the :math:`\alphathm` entry
  in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_hfl` (default :nml_v:`1.`)
  Scaling factor for horizontal flux perturbations (see the :math:`\alphahfl`
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_gam` (default :nml_v:`1.`)
  Scaling factor for g-mode isolation (see the :math:`\alphagam` term in
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_pi` (default :nml_v:`1.`)
  Scaling factor for p-mode isolation (see the :math:`\alphapi` term in
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_kar` (default :nml_v:`1.`)
  Scaling factor for opacity density partial derivative (see the :math:`\alphakar`
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_kat` (default :nml_v:`1.`)
  Scaling factor for opacity temperature partial derivative (see the :math:`\alphakat`
  entry in the :ref:`osc-physics-switches` section)

:nml_n:`alpha_rht` (default :nml_v:`0.`)
  Scaling factor for time-dependent term in radiative heat equation (see the
  :math:`\alpharht` entry in the :ref:`osc-physics-switches` section)
  
:nml_n:`alpha_trb` (default :nml_v:`0.`)
   Scaling factor for the turbulent mixing length (see the
   :math:`\alphatrb` entry in the :ref:`osc-physics-switches`
   section)

:nml_n:`inertia_norm` (default :nml_v:`'BOTH'`)
  Inertia normalization factor; one of

  - :nml_v:`'RADIAL'` : Radial amplitude squared, :math:`|\xi_{\rm r}|^{2}`, evaluated at :nml_v:`x_ref`
  - :nml_v:`'HORIZ'` : Horizontal amplitude squared, :math:`|\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_v:`x_ref`
  - :nml_v:`'BOTH'` : Overall amplitude squared, :math:`|\xi_{\rm r}|^{2} + |\lambda| |\xi_{\rm h}|^{2}`, evaluated at :nml_v:`x_ref`

:nml_n:`time_factor` (default :nml_v:`'OSC'`)
  Time-dependence factor in pulsation equations; one of:

  - :nml_v:`'OSC'` : Oscillatory, :math:`\propto \exp(-{\rm i} \sigma t)`
  - :nml_v:`'EXP'` : Exponential, :math:`\propto \exp(-\sigma t)`

:nml_n:`conv_scheme` (default :nml_v:`'FROZEN_PESNELL_1'`)
  Scheme for treating convection; one of:

  - :nml_v:`'FROZEN_PESNELL_1'` : Freeze convective heating altogether;
    case 1 described by :ads_citet:`pesnell:1990`
  - :nml_v:`'FROZEN_PESNELL_4'` : Freeze Lagrangian perturbation of convective luminosity;
    case 4 described by :ads_citet:`pesnell:1990`

:nml_n:`deps_scheme` (default :nml_v:`'MODEL'`)
  Scheme for calculating nuclear energy generation partials :math:`\epsnucrho` and :math:`\epsnucT`; one of:

  - :nml_v:`'MODEL'` : Use values from model
  - :nml_v:`'FILE'` : Use complex (phase-lagged) values from separate file

:nml_n:`deps_file` (default :nml_v:`''`)
  Name of epsilon partial derivatives file (when :nml_n:`deps_scheme`\ =\ :nml_v:`'FILE'`)

:nml_n:`deps_file_format` (default :nml_v:`'WOLF'`)
  Format of epsilon partial derivative file (when :nml_n:`deps_scheme`\ =\ :nml_v:`'FILE'`); one of:

  - :nml_v:`'WOLF'` : Format used in preparation of :ads_citet:`wolf:2018`

:nml_n:`x_ref` (default :nml_v:`1` or outer grid point, whichever is smaller)
  Reference fractional radius for photosphere, normalizations etc.

:nml_n:`x_atm` (default :nml_v:`-1`, implying outer grid point)
  Fractional radius for convection-zone crossover point of :math:`\pi/\gamma` modes (isolated p and g modes; see :ads_citealp:`ong:2020`)
   
:nml_n:`adiabatic` (default :nml_v:`.TRUE.`)
  Flag to perform adiabatic calculations
  
:nml_n:`nonadiabatic` (default :nml_v:`.FALSE.`)
  Flag to perform non-adiabatic calculations
  
:nml_n:`quasiad_eigfuncs` (default :nml_v:`.FALSE.`)
  Flag to calculate quasi-adiabatic entropy/luminosity eigenfunctions
  during adiabatic calculations

:nml_n:`reduce_order` (default :nml_v:`.TRUE.`)
   Flag to reduce the order of the *adiabatic* radial-pulsation
   equations from 4 to 2

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
