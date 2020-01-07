Mode Files
==========

A mode file gathers together information about a single mode found
during a GYRE run. The data written to a mode file is controlled by
the :nml_o:`mode_item_list` parameter of the :nml_g:`ad_output`
namelist group (for adiabatic calculations) and the
:nml_g:`nad_output` namelist group (for nonadiabatic
calculations). This parameter is a comma-separated list of items to
appear in the mode files; the default value is
:nml_l:`'l,n_pg,omega,freq,x,xi_r,xi_h'`.

The items come in two flavors:

- *scalar* items comprise a single value, typically pertaining either
  to the star as a whole (i.e., a global quantity) or to a specific
  location in the star

- *array* items comprise a sequence of values, with each value
  pertaining to a single point in the discrete grid used to solve the
  oscillation equations. The sequence runs from the inner boundary to
  the outer boundary

The following subsections describe the items that may appear in a
:nml_o:`mode_item_list` parameter, grouped together by functional
area.

Solution Data
-------------

:nml_l:`n` (integer scalar)
  Number of grid points :math:`n`
  
:nml_l:`x` (real array)
  Fractional radius :math:`x \equiv r/R`

:nml_l:`y_1` (complex array)
  Solution variable :math:`y_{1}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`y_2` (complex array)
  Solution variable :math:`y_{2}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`y_3` (complex array)
  Solution variable :math:`y_{3}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`y_4` (complex array)
  Solution variable :math:`y_{4}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`y_5` (complex array)
  Solution variable :math:`y_{5}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`y_6` (complex array)
  Solution variable :math:`y_{6}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`omega` (complex scalar)
  Dimensionless eigenfrequency :math:`\omega`

Observables
-----------

:nml_l:`freq` (complex scalar)
  Dimensioned eigenfrequency. The units and reference frame are
  controlled by :nml_o:`freq_units` and :nml_o:`freq_frame` parameters
  of the :nml_g:`ad_output` and :nml_g:`nad_output` namelist groups
       
:nml_l:`freq_units` (character scalar)
  Units of :nml_l:`freq`

:nml_l:`freq_frame` (character scalar)
  Reference frame of :nml_l:`freq`

:nml_l:`f_T` (real scalar)
  Effective temperature perturbation amplitude :math:`f_{\rm
  T}`. Evaluated using eqn. 5 of :cite:`Dupret:2003a`

:nml_l:`f_g` (real scalar)
  Effective gravity perturbation amplitude :math:`f_{\rm
  g}`. Evaluated using eqn.  6 of :cite:`Dupret:2003a`

:nml_l:`psi_T` (real scalar)
  Effective temperature perturbation phase :math:`\psi_{\rm
  T}`. Evaluated using eqn. 5 of :cite:`Dupret:2003a`

:nml_l:`psi_g` (real scalar)
  Effective gravity perturbation phase :math:`\psi_{\rm g}`

Classification & Validation
---------------------------

:nml_l:`j` (integer scalar)
  Unique mode index :math:`j`. The first mode
  found during the GYRE run has :math:`j=1`, the second
  :math:`j=2`, and so on
  
:nml_l:`l` (integer scalar)
  Harmonic degree :math:`\ell`

:nml_l:`l_i` (complex scalar)
  Effective harmonic degree at inner boundary :math:`\ell_{\rm i}` 

:nml_l:`m` (integer scalar)
  Azimuthal order :math:`m`

:nml_l:`n_p` (integer scalar)
  Acoustic-wave winding number :math:`n_{\rm p}`
  
:nml_l:`n_g` (integer scalar)
  Gravity-wave winding number :math:`n_{\rm g}`

:nml_l:`n_pg` (integer scalar)
  Radial order :math:`n_{\rm pg}` within the Eckart-Scuflaire-Osaki-Takata
  scheme (see :cite:`Takata:2006b`)
   
:nml_l:`omega_int` (complex scalar)
  Dimensionless eigenfrequency :math:`\omega` from integral
  expression. Evaluated using eqn. 1.71 of :cite:`Dupret:2002a`

:nml_l:`Yt_1` (complex array)
  Primary eigenfunction for Takata classification
  :math:`\mathcal{Y}_{1}`. Evaluated using a rescaled eqn. 69 of
  :cite:`Takata:2006b`

:nml_l:`Yt_2` (complex array)
  Secondary eigenfunction for Takata
  classification :math:`\mathcal{Y}_{2}`. Evaluated using a rescaled eqn. 70
  of :cite:`Takata:2006b`

:nml_l:`I_0` (complex array)
  First integral for radial modes :math:`I_{0}`. Evaluated using
  eqn. 42 of :cite:`Takata:2006a`
  
:nml_l:`I_1` (complex array)
  First integral for dipole modes :math:`I_{1}`. Evaluated using
  eqn. 43 of :cite:`Takata:2006a`
  
:nml_l:`prop_type` (complex array)
  Propagation type :math:`\varpi` based on local dispersion
  relation. :math:`\varpi = 1` in acoustic-wave regions,
  :math:`\varpi=-1` in gravity-wave regions, and :math:`\varpi=0` in
  evanescent regions

Perturbations
-------------

:nml_l:`x_ref` (real scalar)
  Fractional radius of reference location :math:`x_{\rm ref}`

:nml_l:`xi_r_ref` (complex scalar)
  Radial displacement perturbation :math:`\xi_{\rm r}` at reference location
  :math:`x_{\rm ref}`, in units of :math:`R`

:nml_l:`xi_h_ref` (complex scalar)
  Horizontal displacement perturbation :math:`\xi_{\rm h}` at reference
  location :math:`x_{\rm ref}`, in units of :math:`R`

:nml_l:`eul_phi_ref` (complex scalar)
  Eulerian potential perturbation :math:`\Phi'` at reference location
  :math:`x_{\rm ref}`, in units of :math:`G M/R`

:nml_l:`deul_phi_ref` (complex scalar)
  Eulerian potential gradient perturbation :math:`{\rm d}\Phi'/{\rm d}x` at
  reference location :math:`x_{\rm ref}`, in units of :math:`G M/R^{2}`

:nml_l:`lag_S_ref` (complex scalar)
  Lagrangian specific entropy perturbation :math:`\delta S` at
  reference location :math:`x_{\rm ref}`, in units of :math:`c_{P}`

:nml_l:`lag_L_ref` (complex scalar)
  Lagrangian radiative luminosity perturbation :math:`\delta L_{r,{\rm
  R}}` at reference location :math:`x_{\rm ref}`, in units of :math:`L`

:nml_l:`xi_r` (complex array)
  Radial displacement perturbation :math:`\xi_{\rm r}`, in units of
  :math:`R`

:nml_l:`xi_h` (complex array)
  Horizontal displacement perturbation :math:`\xi_{\rm h}`, in units
  of :math:`R`

:nml_l:`eul_phi` (complex array)
  Eulerian potential perturbation :math:`\Phi'`, in units of :math:`G
  M/R`

:nml_l:`deul_phi` (complex array)
  Eulerian potential gradient perturbation :math:`{\rm d}\Phi'/{\rm
  d}x`, in units of :math:`G M/R^{2}`

:nml_l:`lag_S` (complex array)
  Lagrangian specific entropy perturbation :math:`\delta S`, in units
  of :math:`c_{P}`

:nml_l:`lag_L` (complex array)
  Lagrangian radiative luminosity peturbation :math:`\delta L_{r,{\rm
  R}}`, in units of :math:`L`

:nml_l:`eul_P` (complex array)
  Eulerian total pressure perturbation :math:`P'`, in units of
  :math:`P`

:nml_l:`eul_rho` (complex array)
  Eulerian density perturbation :math:`\rho'`, in units of
  :math:`\rho`

:nml_l:`eul_T` (complex array)
  Eulerian temperature perturbation :math:`T'`, in units of :math:`T`
       
:nml_l:`lag_P` (complex array)
  Lagrangian total pressure perturbation :math:`\delta P`, in units of
  :math:`P`

:nml_l:`lag_rho` (complex array)
  Lagrangian density perturbation :math:`\delta \rho`, in units of
  :math:`\rho`

:nml_l:`lag_T` (complex array)
  Lagrangian temperature perturbation :math:`\delta T`, in units of
  :math:`T`

Energetics & Transport
----------------------

:nml_l:`eta` (real scalar)
  Normalized growth rate :math:`\eta`. Evaluated using expression in
  text of page 1186 of :cite:`Stellingwerf:1978a`

:nml_l:`E` : (real scalar)
  Mode inertia :math:`E`, in units of :math:`M R^{2}`. Evaluated
  by integrating :math:`{\rm d}E/{\rm d}x`

:nml_l:`E_p` (real scalar)
  Acoustic inertia :math:`E_{\rm p}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  acoustic-wave propagation regions

:nml_l:`E_g` (real scalar)
  Gravity inertia :math:`E_{\rm g}`, in units of :math:`M
  R^{2}`. Evaluated by integrating :math:`{\rm d}E/{\rm d}x` in
  gravity-wave propagation regions

:nml_l:`E_norm` (real scalar)
  Normalized inertia :math:`E_{\rm norm}`. The normalization is
  controlled by the :nml_o:`inertia_norm` parameter of the
  :nml_g:`osc` namelist group

:nml_l:`E_ratio` (real scalar)
  Ratio of mode inertia inside/outside the reference location
  :math:`x_{\rm ref}`

:nml_l:`H` (real scalar)
  Mode energy :math:`H`, in units of :math:`G M^{2}/R`

:nml_l:`W` (real scalar)
  Mode work :math:`W`, in units of :math:`G M^{2}/R`. Evaluated by
  integrating :math:`{\rm d}W/{\rm d}x`

:nml_l:`W_eps` (real scalar)
  Mode nuclear work :math:`W_{\epsilon}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}W_{\epsilon}/{\rm
  d}x`

:nml_l:`tau_ss` (real scalar)
  Steady-state mode torque :math:`\tau_{\rm ss}`, in units of :math:`G
  M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm ss}/{\rm
  d}x`

:nml_l:`tau_tr` (real scalar)
  Transient total mode torque :math:`\tau_{\rm tr}`, in units of
  :math:`G M^{2}/R`. Evaluated by integrating :math:`{\rm d}\tau_{\rm
  tr}/{\rm d}x`

:nml_l:`dE_dx` (real array)
  Differential inertia :math:`{\rm d}E/{\rm d}x`, in units of :math:`M
  R^{2}`

:nml_l:`dW_dx` (real array)
  Differential work :math:`{\rm d}W/{\rm d}x`, in units of :math:`G
  M^{2}/R`. Evaluated using eqn. 25.9 of :cite:`Unno:1989a`

:nml_l:`dW_eps_dx` (real array)
  Differential nuclear work :math:`{\rm d}W_{epsilon}/{\rm d}x`,
  in units of :math:`G M^{2}/R`. Evaluated using eqn. 25.9 of
  :cite:`Unno:1989a`

:nml_l:`dtau_dx_ss` (real array)
  Steady-state differential torque :math:`{\rm d}\tau_{\rm ss}/{\rm
  d}x`, in units of :math:`G M^{2}/R`

:nml_l:`dtau_dx_tr` (real array)
  Transient differential torque :math:`{\rm d}\tau_{\rm tr}/{\rm d}x`,
  in units of :math:`G M^{2}/R`

:nml_l:`alpha_0` (real array)
  Excitation coefficient :math:`\alpha_{0}`. Evaluated using eqn. 26.10
  of :cite:`Unno:1989a`

:nml_l:`alpha_1` (real array)
  Excitation coefficient :math:`\alpha_{1}`. Evaluated using eqn. 26.12
  of :cite:`Unno:1989a`

Rotation
--------

:nml_l:`beta` (real scalar)
  Rotation splitting coefficient :math:`\beta`. Evaluated by
  integrating :math:`{\rm d}\beta/{\rm d}x`

:nml_l:`dbeta_dx` (real array)
  Unnormalized rotation splitting kernel :math:`{\rm d}\beta/{\rm
  d}x`. Evaluated using eqn. 3.357 of :cite:`Aerts:2010a`
	 
:nml_l:`lambda` (complex array)
  Angular eigenvalue :math:`\lambda`

Stellar Structure
-----------------

:nml_l:`M_star` (real scalar)
  stellar mass, in units of :math:`{\rm g}` [#only_evol]_

:nml_l:`R_star` (real scalar)
  stellar radius, in units of :math:`{\rm cm}` [#only_evol]_

:nml_l:`L_star` (real scalar)
  stellar luminosity, in units of :math:`{\rm erg\,s^{-1}}` [#only_evol]_

:nml_l:`Delta_p` (real scalar)
  Asymptotic p-mode large frequency separation :math:`\Delta \nu`,
  in units of :math:`\sqrt{GM/R^{3}}`

:nml_l:`Delta_g` (real scalar)
  Asymptotic g-mode inverse period separation :math:`(\Delta
  P)^{-1}`, in units of :math:`\sqrt{GM/R^{3}}`

:nml_l:`V_2` (real array)
  Dimensionless structure coefficient :math:`V_{2}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`As` (real array)
  Dimensionless structure coefficient :math:`A^{*}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_l:`U` (real array)
  Dimensionless structure coefficient :math:`U`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
   
:nml_l:`c_1` (real array)
  Dimensionless structure coefficient :math:`c_{1}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`Gamma_1` (real array)
  Adiabatic exponent :math:`\Gamma_{1}`,
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`nabla` (real array)
  Dimensionless temperature gradient :math:`\nabla`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_l:`nabla_ad` (real array)
  Adiabatic tempertature gradient :math:`\nabla_{\rm ad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`dnabla_ad` (real array)
  Dimensionless gradient :math:`\partial \nabla_{\rm ad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
   
:nml_l:`\delta` (real array)
  Thermodynamic coefficient :math:`\delta`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`c_lum` (real array)
  Dimensionless structure coefficient :math:`c_{\rm lum}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`c_rad` (real array)
  Dimensionless structure coefficient :math:`c_{\rm rad}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`c_thn` (real array)
  Dimensionless structure coefficient :math:`c_{\rm thn}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`c_thk` (real array)
  Dimensionless structure coefficient :math:`c_{\rm thk}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`c_eps` (real array)
  Dimensionless structure coefficient :math:`c_{\epsilon}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`eps_rho` (real array)
  Energy generation partial :math:`\epsilon_{\rho}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`
  
:nml_l:`eps_T` (real array)
  Energy generation partial :math:`\epsilon_{T}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`kap_rho` (real array)
  Opacity partial :math:`\kappa_{\rho}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`kap_T` (real array)
  Opacity partial :math:`\kappa_{T}`, defined in
  :repo:`equations.pdf <doc/equations.pdf>`

:nml_l:`Omega_rot` (real array)
  Rotation angular frequency, in units of :math:`\sqrt{GM/R^{3}}`

:nml_l:`M_r` (real array)
  Mass coordinate, in units of :math:`{\rm g}` [#only_evol]_

:nml_l:`P` (real array)
  Total pressure, in units of :math:`{\rm dyn\,cm^{-2}}` [#only_evol]_

:nml_l:`\rho` (real array)
  Density, in units of :math:`{\rm g\,cm^{-3}}` [#only_evol]_

:nml_l:`T` (real array)
  Temperature, in units of :math:`{\rm K}` [#only_evol]_

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_o:`model_type`
                is :nml_l:`'EVOL'`
