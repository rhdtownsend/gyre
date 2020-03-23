.. _mode-files:

Mode Files
==========

A mode file gathers together information about a single mode found
during a GYRE run. The data written to a mode file is controlled by
the :nml_n:`mode_item_list` parameter of the :nml_g:`ad_output`
namelist group (for adiabatic calculations) and the
:nml_g:`nad_output` namelist group (for nonadiabatic
calculations). This parameter is a comma-separated list of items to
appear in the mode files. The items come in two flavors:

* *scalar* items comprise a single value, typically pertaining either
  to the star as a whole (i.e., a global quantity) or to a specific
  location in the star

* *array* items comprise a sequence of values, with each value
  pertaining to a single point in the discrete grid used to solve the
  oscillation equations. The sequence runs from the inner boundary to
  the outer boundary

The following subsections describe the items that may appear in a
:nml_n:`mode_item_list` parameter, grouped together by functional
area.

Solution Data
-------------

:nml_v:`n` (integer scalar)
  Number of grid points :math:`n`
  
:nml_v:`x` (real array)
  Fractional radius :math:`x \equiv r/R`

:nml_v:`y_1` (complex array)
  Solution variable :math:`y_{1}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`y_2` (complex array)
  Solution variable :math:`y_{2}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`y_3` (complex array)
  Solution variable :math:`y_{3}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`y_4` (complex array)
  Solution variable :math:`y_{4}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`y_5` (complex array)
  Solution variable :math:`y_{5}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`y_6` (complex array)
  Solution variable :math:`y_{6}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`omega` (complex scalar)
  Dimensionless eigenfrequency :math:`\omega`

Observables
-----------

:nml_v:`freq` (complex scalar)
  Dimensioned eigenfrequency. The units and reference frame are
  controlled by :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
  of the :nml_g:`ad_output` and :nml_g:`nad_output` namelist groups
       
:nml_v:`freq_units` (character scalar)
  Units of :nml_v:`freq`

:nml_v:`freq_frame` (character scalar)
  Reference frame of :nml_v:`freq`

:nml_v:`f_T` (real scalar)
  Effective temperature perturbation amplitude :math:`f_{\rm
  T}`. Evaluated using eqn. 5 of :ads_citet:`dupret:2003`

:nml_v:`f_g` (real scalar)
  Effective gravity perturbation amplitude :math:`f_{\rm
  g}`. Evaluated using eqn.  6 of :ads_citet:`dupret:2003`

:nml_v:`psi_T` (real scalar)
  Effective temperature perturbation phase :math:`\psi_{\rm
  T}`. Evaluated using eqn. 5 of :ads_citet:`dupret:2003`

:nml_v:`psi_g` (real scalar)
  Effective gravity perturbation phase :math:`\psi_{\rm g}`

Classification & Validation
---------------------------

:nml_v:`j` (integer scalar)
  Unique mode index :math:`j`. The first mode
  found during the GYRE run has :math:`j=1`, the second
  :math:`j=2`, and so on
  
:nml_v:`l` (integer scalar)
  Harmonic degree :math:`\ell`

:nml_v:`l_i` (complex scalar)
  Effective harmonic degree at inner boundary :math:`\ell_{\rm i}` 

:nml_v:`m` (integer scalar)
  Azimuthal order :math:`m`

:nml_v:`n_p` (integer scalar)
  Acoustic-wave winding number :math:`n_{\rm p}`
  
:nml_v:`n_g` (integer scalar)
  Gravity-wave winding number :math:`n_{\rm g}`

:nml_v:`n_pg` (integer scalar)
  Radial order :math:`n_{\rm pg}` within the Eckart-Scuflaire-Osaki-Takata
  scheme (see :ads_citealp:`takata:2006b`)
   
:nml_v:`omega_int` (complex scalar)
  Dimensionless eigenfrequency :math:`\omega` from integral
  expression. Evaluated using eqn. 1.71 of Marc-Antoine's Dupret's PhD thesis

:nml_v:`Yt_1` (complex array)
  Primary eigenfunction for Takata classification
  :math:`\mathcal{Y}_{1}`. Evaluated using a rescaled eqn. 69 of
  :ads_citet:`takata:2006b`

:nml_v:`Yt_2` (complex array)
  Secondary eigenfunction for Takata
  classification :math:`\mathcal{Y}_{2}`. Evaluated using a rescaled eqn. 70
  of :ads_citet:`takata:2006b`

:nml_v:`I_0` (complex array)
  First integral for radial modes :math:`I_{0}`. Evaluated using
  eqn. 42 of :ads_citet:`takata:2006a`
  
:nml_v:`I_1` (complex array)
  First integral for dipole modes :math:`I_{1}`. Evaluated using
  eqn. 43 of :ads_citet:`takata:2006a`
  
:nml_v:`prop_type` (complex array)
  Propagation type :math:`\varpi` based on local dispersion
  relation. :math:`\varpi = 1` in acoustic-wave regions,
  :math:`\varpi=-1` in gravity-wave regions, and :math:`\varpi=0` in
  evanescent regions

Perturbations
-------------

:nml_v:`x_ref` (real scalar)
  Fractional radius of reference location :math:`x_{\rm ref}`

:nml_v:`xi_r_ref` (complex scalar)
  Radial displacement perturbation :math:`\xi_{\rm r}` at reference location
  :math:`x_{\rm ref}`, in units of :math:`R`

:nml_v:`xi_h_ref` (complex scalar)
  Horizontal displacement perturbation :math:`\xi_{\rm h}` at reference
  location :math:`x_{\rm ref}`, in units of :math:`R`

:nml_v:`eul_phi_ref` (complex scalar)
  Eulerian potential perturbation :math:`\Phi'` at reference location
  :math:`x_{\rm ref}`, in units of :math:`G M/R`

:nml_v:`deul_phi_ref` (complex scalar)
  Eulerian potential gradient perturbation :math:`{\rm d}\Phi'/{\rm d}x` at
  reference location :math:`x_{\rm ref}`, in units of :math:`G M/R^{2}`

:nml_v:`lag_S_ref` (complex scalar)
  Lagrangian specific entropy perturbation :math:`\delta S` at
  reference location :math:`x_{\rm ref}`, in units of :math:`c_{P}`

:nml_v:`lag_L_ref` (complex scalar)
  Lagrangian radiative luminosity perturbation :math:`\delta L_{r,{\rm
  R}}` at reference location :math:`x_{\rm ref}`, in units of :math:`L`

:nml_v:`xi_r` (complex array)
  Radial displacement perturbation :math:`\xi_{\rm r}`, in units of
  :math:`R`

:nml_v:`xi_h` (complex array)
  Horizontal displacement perturbation :math:`\xi_{\rm h}`, in units
  of :math:`R`

:nml_v:`eul_phi` (complex array)
  Eulerian potential perturbation :math:`\Phi'`, in units of :math:`G
  M/R`

:nml_v:`deul_phi` (complex array)
  Eulerian potential gradient perturbation :math:`{\rm d}\Phi'/{\rm
  d}x`, in units of :math:`G M/R^{2}`

:nml_v:`lag_S` (complex array)
  Lagrangian specific entropy perturbation :math:`\delta S`, in units
  of :math:`c_{P}`

:nml_v:`lag_L` (complex array)
  Lagrangian radiative luminosity peturbation :math:`\delta L_{r,{\rm
  R}}`, in units of :math:`L`

:nml_v:`eul_P` (complex array)
  Eulerian total pressure perturbation :math:`P'`, in units of
  :math:`P`

:nml_v:`eul_rho` (complex array)
  Eulerian density perturbation :math:`\rho'`, in units of
  :math:`\rho`

:nml_v:`eul_T` (complex array)
  Eulerian temperature perturbation :math:`T'`, in units of :math:`T`
       
:nml_v:`lag_P` (complex array)
  Lagrangian total pressure perturbation :math:`\delta P`, in units of
  :math:`P`

:nml_v:`lag_rho` (complex array)
  Lagrangian density perturbation :math:`\delta \rho`, in units of
  :math:`\rho`

:nml_v:`lag_T` (complex array)
  Lagrangian temperature perturbation :math:`\delta T`, in units of
  :math:`T`

Energetics & Transport
----------------------

:nml_v:`eta` (real scalar)
  Normalized growth rate :math:`\eta`. Evaluated using expression in
  text of page 1186 of :ads_citet:`stellingwerf:1978`

:nml_v:`E` : (real scalar)
  Mode inertia :math:`E`, in units of :math:`M R^{2}`. Evaluated
  by integrating :math:`{\rm d}E/{\rm d}x`

:nml_v:`E_p` (real scalar)
  Acoustic inertia :math:`E_{\rm p}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  acoustic-wave propagation regions

:nml_v:`E_g` (real scalar)
  Gravity inertia :math:`E_{\rm g}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  gravity-wave propagation regions

:nml_v:`E_norm` (real scalar)
  Normalized inertia :math:`E_{\rm norm}`. The normalization is
  controlled by the :nml_n:`inertia_norm` parameter of the
  :nml_g:`osc` namelist group

:nml_v:`E_ratio` (real scalar)
  Ratio of mode inertia inside/outside the reference location
  :math:`x_{\rm ref}`

:nml_v:`H` (real scalar)
  Mode energy :math:`H`, in units of :math:`G M^{2}/R`

:nml_v:`W` (real scalar)
  Mode work :math:`W`, in units of :math:`G M^{2}/R`. Evaluated by
  integrating :math:`{\rm d}W/{\rm d}x`

:nml_v:`W_eps` (real scalar)
  Mode nuclear work :math:`W_{\epsilon}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}W_{\epsilon}/{\rm
  d}x`

:nml_v:`tau_ss` (real scalar)
  Steady-state mode torque :math:`\tau_{\rm ss}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm ss}/{\rm
  d}x`

:nml_v:`tau_tr` (real scalar)
  Transient total mode torque :math:`\tau_{\rm tr}`, in units of
  :math:`G M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm
  tr}/{\rm d}x`

:nml_v:`dE_dx` (real array)
  Differential inertia :math:`{\rm d}E/{\rm d}x`, in units of :math:`M
  R^{2}`

:nml_v:`dW_dx` (real array)
  Differential work :math:`{\rm d}W/{\rm d}x`, in units of :math:`G
  M^{2}/R`. Evaluated using eqn. 25.9 of :ads_citet:`unno:1989`

:nml_v:`dW_eps_dx` (real array)
  Differential nuclear work :math:`{\rm d}W_{epsilon}/{\rm d}x`,
  in units of :math:`G M^{2}/R`. Evaluated using eqn. 25.9 of
  :ads_citet:`unno:1989`

:nml_v:`dtau_dx_ss` (real array)
  Steady-state differential torque :math:`{\rm d}\tau_{\rm ss}/{\rm
  d}x`, in units of :math:`G M^{2}/R`

:nml_v:`dtau_dx_tr` (real array)
  Transient differential torque :math:`{\rm d}\tau_{\rm tr}/{\rm d}x`,
  in units of :math:`G M^{2}/R`

:nml_v:`alpha_0` (real array)
  Excitation coefficient :math:`\alpha_{0}`. Evaluated using eqn. 26.10
  of :ads_citet:`unno:1989`

:nml_v:`alpha_1` (real array)
  Excitation coefficient :math:`\alpha_{1}`. Evaluated using eqn. 26.12
  of :ads_citet:`unno:1989`

Rotation
--------

:nml_v:`beta` (real scalar)
  Rotation splitting coefficient :math:`\beta`. Evaluated by
  integrating :math:`{\rm d}\beta/{\rm d}x`

:nml_v:`dbeta_dx` (real array)
  Unnormalized rotation splitting kernel :math:`{\rm d}\beta/{\rm
  d}x`. Evaluated using eqn. 3.357 of :ads_citet:`aerts:2010`
	 
:nml_v:`lambda` (complex array)
  Angular eigenvalue :math:`\lambda`

Stellar Structure
-----------------

:nml_v:`M_star` (real scalar)
  stellar mass, in units of :math:`{\rm g}` [#only_evol]_

:nml_v:`R_star` (real scalar)
  stellar radius, in units of :math:`{\rm cm}` [#only_evol]_

:nml_v:`L_star` (real scalar)
  stellar luminosity, in units of :math:`{\rm erg\,s^{-1}}` [#only_evol]_

:nml_v:`Delta_p` (real scalar)
  Asymptotic p-mode large frequency separation :math:`\Delta \nu`,
  in units of :math:`\sqrt{GM/R^{3}}`

:nml_v:`Delta_g` (real scalar)
  Asymptotic g-mode inverse period separation :math:`(\Delta
  P)^{-1}`, in units of :math:`\sqrt{GM/R^{3}}`

:nml_v:`V_2` (real array)
  Dimensionless structure coefficient :math:`V_{2}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`As` (real array)
  Dimensionless structure coefficient :math:`A^{*}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_v:`U` (real array)
  Dimensionless structure coefficient :math:`U`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
   
:nml_v:`c_1` (real array)
  Dimensionless structure coefficient :math:`c_{1}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`Gamma_1` (real array)
  Adiabatic exponent :math:`\Gamma_{1}`,
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`nabla` (real array)
  Dimensionless temperature gradient :math:`\nabla`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_v:`nabla_ad` (real array)
  Adiabatic tempertature gradient :math:`\nabla_{\rm ad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`dnabla_ad` (real array)
  Dimensionless gradient :math:`\partial \nabla_{\rm ad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
   
:nml_v:`\delta` (real array)
  Thermodynamic coefficient :math:`\delta`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`c_lum` (real array)
  Dimensionless structure coefficient :math:`c_{\rm lum}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`c_rad` (real array)
  Dimensionless structure coefficient :math:`c_{\rm rad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`c_thn` (real array)
  Dimensionless structure coefficient :math:`c_{\rm thn}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`c_thk` (real array)
  Dimensionless structure coefficient :math:`c_{\rm thk}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`c_eps` (real array)
  Dimensionless structure coefficient :math:`c_{\epsilon}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`eps_rho` (real array)
  Energy generation partial :math:`\epsilon_{\rho}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_v:`eps_T` (real array)
  Energy generation partial :math:`\epsilon_{T}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`kap_rho` (real array)
  Opacity partial :math:`\kappa_{\rho}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`kap_T` (real array)
  Opacity partial :math:`\kappa_{T}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_v:`Omega_rot` (real array)
  Rotation angular frequency, in units of :math:`\sqrt{GM/R^{3}}`

:nml_v:`M_r` (real array)
  Mass coordinate, in units of :math:`{\rm g}` [#only_evol]_

:nml_v:`P` (real array)
  Total pressure, in units of :math:`{\rm dyn\,cm^{-2}}` [#only_evol]_

:nml_v:`\rho` (real array)
  Density, in units of :math:`{\rm g\,cm^{-3}}` [#only_evol]_

:nml_v:`T` (real array)
  Temperature, in units of :math:`{\rm K}` [#only_evol]_

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_n:`model_type`
                is :nml_v:`'EVOL'`
