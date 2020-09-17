.. _summary-files:

Summary Files
=============

A summary file gathers together information about all modes found during
a GYRE run. The data written to summary files is controlled by the
:nml_n:`summary_item_list` parameter of the :nml_g:`ad_output`
namelist group (for adiabatic calculations) and the
:nml_g:`nad_output` namelist group (for nonadiabatic
calculations). This parameter is a comma-separated list of items to
appear in the summary file. The items come in two flavors:

* *scalar* items comprise a single value, typically pertaining to all
  modes (i.e., a global quantity)

* *array* items comprise a sequence of values, with each value
  pertaining to a single mode. The sequence follows the same order in
  which modes were found during the GYRE run.

The following subsections describe the items that may appear in a
:nml_n:`summary_item_list` parameter, grouped together by functional
area.

Solution Data
-------------

:nml_v:`omega` (complex array)
  Dimensionless eigenfrequency :math:`\omega`

Observables
-----------

:nml_v:`freq` (complex array)
  Dimensioned eigenfrequency. The units and reference frame are
  controlled by :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
  of the :nml_g:`ad_output` and :nml_g:`nad_output` namelist groups
       
:nml_v:`freq_units` (character scalar)
  Units of :nml_v:`freq`

:nml_v:`freq_frame` (character scalar)
  Reference frame of :nml_v:`freq`

:nml_v:`f_T` (real array)
  Effective temperature perturbation amplitude :math:`f_{\rm
  T}`. Evaluated using eqn. 5 of :ads_citet:`dupret:2003`

:nml_v:`f_g` (real array)
  Effective gravity perturbation amplitude :math:`f_{\rm
  g}`. Evaluated using eqn.  6 of :ads_citet:`dupret:2003`

:nml_v:`psi_T` (real array)
  Effective temperature perturbation phase :math:`\psi_{\rm
  T}`. Evaluated using eqn. 5 of :ads_citet:`dupret:2003`

:nml_v:`psi_g` (real array)
  Effective gravity perturbation phase :math:`\psi_{\rm g}`

Classification & Validation
---------------------------

:nml_v:`j` (integer array)
  Unique mode index :math:`j`. The first mode
  found during the GYRE run has :math:`j=1`, the second
  :math:`j=2`, and so on
  
:nml_v:`l` (integer array)
  Harmonic degree :math:`\ell`
  
:nml_v:`l_i` (complex array)
  Effective harmonic degree at inner boundary :math:`\ell_{\rm i}` 

:nml_v:`m` (integer array)
  Azimuthal order :math:`m`
  
:nml_v:`n_p` (integer array)
  Acoustic-wave winding number :math:`n_{\rm p}`
  
:nml_v:`n_g` (integer array)
  Gravity-wave winding number :math:`n_{\rm g}`

:nml_v:`n_pg` (integer array)
  Radial order :math:`n_{\rm pg}` within the Eckart-Scuflaire-Osaki-Takata
  scheme (see :ads_citealp:`takata:2006b`)
   
:nml_v:`omega_int` (complex array)
  Dimensionless eigenfrequency :math:`\omega` from integral
  expression. Evaluated using eqn. 1.71 of Marc-Antoine Dupret's PhD thesis

Perturbations
-------------
  
:nml_v:`x_ref` (real array)
  Fractional radius of reference location :math:`x_{\rm ref}`

:nml_v:`xi_r_ref` (complex array)
  Radial displacement perturbation :math:`\xi_{\rm r}` at reference location
  :math:`x_{\rm ref}`, in units of :math:`R`

:nml_v:`xi_h_ref` (complex array)
  Horizontal displacement perturbation :math:`\xi_{\rm h}` at reference
  location :math:`x_{\rm ref}`, in units of :math:`R`

:nml_v:`eul_phi_ref` (complex array)
  Eulerian potential perturbation :math:`\Phi'` at reference location
  :math:`x_{\rm ref}`, in units of :math:`G M/R`

:nml_v:`deul_phi_ref` (complex array)
  Eulerian potential gradient perturbation :math:`{\rm d}\Phi'/{\rm d}x` at
  reference location :math:`x_{\rm ref}`, in units of :math:`G M/R^{2}`

:nml_v:`lag_S_ref` (complex array)
  Lagrangian specific entropy perturbation :math:`\delta S` at
  reference location :math:`x_{\rm ref}`, in units of :math:`c_{P}`

:nml_v:`lag_L_ref` (complex array)
  Lagrangian radiative luminosity perturbation :math:`\delta L_{r,{\rm
  R}}` at reference location :math:`x_{\rm ref}`, in units of :math:`L`

Energetics & Transport
----------------------

:nml_v:`eta` (real array)
  Normalized growth rate :math:`\eta`. Evaluated using expression in
  text of page 1186 of :ads_citet:`stellingwerf:1978`

:nml_v:`E` : (real array)
  Mode inertia :math:`E`, in units of :math:`M R^{2}`. Evaluated
  by integrating :math:`{\rm d}E/{\rm d}x`

:nml_v:`E_p` (real array)
  Acoustic inertia :math:`E_{\rm p}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  acoustic-wave propagation regions

:nml_v:`E_p` (real array)
  Gravity inertia :math:`E_{\rm g}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  gravity-wave propagation regions

:nml_v:`E_norm` (real array)
  Normalized inertia :math:`E_{\rm norm}`. The normalization is
  controlled by the :nml_n:`inertia_norm` parameter of the
  :nml_g:`osc` namelist group

:nml_v:`E_ratio` (real array)
  Ratio of mode inertia inside/outside the reference location
  :math:`x_{\rm ref}`

:nml_v:`H` (real array)
  Mode energy :math:`H`, in units of :math:`G M^{2}/R`

:nml_v:`W` (real array)
  Mode work :math:`W`, in units of :math:`G M^{2}/R`. Evaluated by
  integrating :math:`{\rm d}W/{\rm d}x`

:nml_v:`W_eps` (real array)
  Mode nuclear work :math:`W_{\epsilon}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}W_{\epsilon}/{\rm
  d}x`

:nml_v:`tau_ss` (real array)
  Steady-state mode torque :math:`\tau_{\rm ss}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm ss}/{\rm
  d}x`

:nml_v:`tau_tr` (real array)
  Transient total mode torque :math:`\tau_{\rm tr}`, in units of
  :math:`G M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm
  tr}/{\rm d}x`

Rotation
--------

:nml_v:`beta` (real array)
  Rotation splitting coefficient :math:`\beta`. Evaluated by
  integrating :math:`{\rm d}\beta/{\rm d}x`

Stellar Structure
-----------------

:nml_v:`M_star` (real scalar)
  stellar mass, in units of :math:`{\rm g}`\ [#only_evol]_

:nml_v:`R_star` (real scalar)
  stellar radius, in units of :math:`{\rm cm}`\ [#only_evol]_

:nml_v:`L_star` (real scalar)
  stellar luminosity, in units of :math:`{\rm erg\,s^{-1}}`\ [#only_evol]_

:nml_v:`Delta_p` (real array)
  Asymptotic p-mode large frequency separation :math:`\Delta \nu`,
  in units of :math:`\sqrt{GM/R^{3}}`

:nml_v:`Delta_g` (real array)
  Asymptotic g-mode inverse period separation :math:`(\Delta
  P)^{-1}`, in units of :math:`\sqrt{GM/R^{3}}`

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_n:`model_type`
                is :nml_v:`'EVOL'`
