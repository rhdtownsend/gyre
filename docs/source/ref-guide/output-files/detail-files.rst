.. _detail-files:

.. ofile:filetype:: detail

.. nml:group:: output
   :no-target:

Detail Files
============

Detail files store spatial quantities, such as eigenfunctions and
differential inertias, for an individual solution (mode, tidal
response, etc.) found during a run. The specific data written to
detail files are controlled by the :nml:option:`detail_item_list`
options of the :nml:group:`ad_output` and :nml:group:`nad_output` namelist
groups (:program:`gyre` adiabatic and non-adiabatic calculations,
respectively) and the :nml:group:`tide_output` namelist group
(:program:`gyre_tides` calculations). These options specify the
data to be written via a comma-separated list of fields.

The following subsections describe the fields that may appear in
:nml:option:`detail_item_list`, grouped together by functional area.

Solution Data
-------------

.. ofile:field:: n
   :type: integer

   Number of spatial grid points :math:`N`

.. ofile:field:: omega
   :type: complex

   Dimensionless angular frequency :math:`\omega`

.. ofile:field:: x
   :type: real
   :dim: :ofile:field:`n`

   Dimensionless radial coordinate :math:`x \equiv r/\Rstar`

.. ofile:field:: x_ref
   :type: real

   Dimensionless radial coordinate :math:`\xref` of reference location

.. ofile:field:: dx_min
   :type: real

   Minimum radial spacing :math:`\Delta x_{\rm min}` in spatial grid

.. ofile:field:: dx_max
   :type: real

   Maximum radial spacing :math:`\Delta x_{\rm max}` in spatial grid

.. ofile:field:: dx_rms
   :type: real

   Root-mean-square spacing :math:`\Delta x_{\rm rms}` of spatial grid

.. ofile:field:: y_1
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{1}`

.. ofile:field:: y_2
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{2}`

.. ofile:field:: y_3
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{3}`

.. ofile:field:: y_4
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{4}`

.. ofile:field:: y_5
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{5}`

.. ofile:field:: y_6
   :type: complex
   :dim: :ofile:field:`n`

   Dependent variable :math:`y_{6}`

.. ofile:field:: chi
   :type: real

   Root-finding convergence parameter :math:`\chi`

.. ofile:field:: n_iter
   :type: integer

   Root-finding number of iterations

The definitions of the dependent variables
:math:`\{y_{1},\ldots,y_{6}\}` are provided in the :ref:`osc-eqns`
chapter.

Observables
-----------

.. ofile:field:: freq
   :type: complex
   :units: controlled by :nml:option:`freq_units` and :nml:option:`freq_frame` options

   Dimensioned frequency

.. ofile:field:: freq_units
   :type: string

   Value of :nml:option:`freq_units` option

.. ofile:field:: freq_frame
   :type: string

   Value of :nml:option:`freq_frame` option

.. ofile:field:: f_T
   :type: real

   Effective temperature perturbation amplitude :math:`f_{T}`;
   evaluated at reference location using eqn. (5) of
   :ads_citet:`dupret:2003`

.. ofile:field:: f_g
   :type: real

   Effective gravity perturbation amplitude :math:`f_{g}`; evaluated
   at reference location using eqn. (6) of :ads_citet:`dupret:2003`

.. ofile:field:: psi_T
   :type: real

   Effective temperature perturbation phase :math:`\psi_{T}`;
   evaluated at reference location using eqn. (5) of
   :ads_citet:`dupret:2003`

.. ofile:field:: psi_g
   :type: real

   Effective gravity perturbation phase :math:`\psi_{g}`; evaluated at
   reference location using eqn. (6) of :ads_citet:`dupret:2003`

Classification & Validation
---------------------------

.. ofile:field:: id
   :type: integer

   Unique mode index

.. ofile:field:: l
   :type: integer

   Harmonic degree :math:`\ell`

.. ofile:field:: l_i
   :type: complex

   Effective harmonic degree at inner boundary :math:`\ell_{\rm i}`

.. ofile:field:: m
   :type: integer

   Azimuthal order :math:`m`

.. ofile:field:: n_p
   :type: integer

   Acoustic-wave winding number :math:`\nump`

.. ofile:field:: n_g
   :type: integer

   Gravity-wave winding number :math:`\numg`

.. ofile:field:: n_pg
   :type: integer

   Radial order :math:`\numpg` within the
   Eckart-Scuflaire-Osaki-Takata scheme (see
   :ads_citealp:`takata:2006b`)

.. ofile:field:: omega_int
   :type: complex

   Dimensionless eigenfrequency :math:`\omega_{\rm int}` based on
   integral expression; evaluated using eqn. (A8) of
   :ads_citet:`townsend:2025`

.. ofile:field:: dzeta_dx
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Frequency weight function :math:`\sderiv{\zeta}{x}`; evaluated from
   the integrand in eqn. (A5) of :ads_citet:`townsend:2025` with
   :math:`n'=n`

.. ofile:field:: zeta
   :type: complex

   Integrated frequency weight :math:`\zeta \equiv \int \sderiv{\zeta}{x} \, \diff{x}`

.. ofile:field:: Yt_1
   :type: complex
   :dim: :ofile:field:`n`

   Primary eigenfunction :math:`\mathcal{Y}_{1}` for Takata
   classification; evaluated using a rescaled eqn. (69) of
   :ads_citet:`takata:2006b`

.. ofile:field:: Yt_2
   :type: complex
   :dim: :ofile:field:`n`

   Secondary eigenfunction :math:`\mathcal{Y}_{2}` for Takata
   classification; evaluated using a rescaled eqn. (70) of
   :ads_citet:`takata:2006b`

.. ofile:field:: I_0
   :type: complex
   :dim: :ofile:field:`n`

   First integral :math:`I_{0}` for radial modes; evaluated using
   eqn. (42) of :ads_citet:`takata:2006a`

.. ofile:field:: I_1
   :type: complex
   :dim: :ofile:field:`n`

   First integral :math:`I_{1}` for dipole modes; evaluated using
   eqn. (43) of :ads_citet:`takata:2006a`

.. ofile:field:: prop_type
   :type: integer
   :dim: :ofile:field:`n`

   Wave propagation type: :math:`\varpi = 1` in acoustic-wave regions,
   :math:`\varpi=-1` in gravity-wave regions, and :math:`\varpi=0` in
   evanescent regions

Perturbations
-------------

.. ofile:field:: xi_r_ref
   :type: complex
   :units: :math:`\Rstar`

   Radial displacement perturbation :math:`\txi_{r,{\rm ref}}` at
   reference location

.. ofile:field:: xi_h_ref
   :type: complex
   :units: :math:`\Rstar`

   Horizontal displacement perturbation :math:`\txi_{\rm h,ref}` at
   reference location

.. ofile:field:: eul_Phi_ref
   :type: complex
   :units:  :math:`G\Mstar/\Rstar`

   Eulerian potential perturbation :math:`\tPhi'_{\rm ref}` at
   reference location

.. ofile:field:: deul_Phi_ref
   :type: complex
   :units: :math:`G\Mstar/\Rstar^{2}`

   Eulerian potential gradient perturbation
   :math:`(\sderiv{\tPhi'}{x})_{\rm ref}` at reference location

.. ofile:field:: lag_S_ref
   :type: complex
   :units: :math:`\cP`

   Lagrangian specific entropy perturbation :math:`\delta\tS_{\rm
   ref}` at reference location

.. ofile:field:: lag_L_ref
   :type: complex
   :units: :math:`\Lstar`

   Lagrangian radiative luminosity perturbation :math:`\delta\tL_{\rm
   R,ref}` at reference location

.. ofile:field:: xi_r
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\Rstar`

   Radial displacement perturbation :math:`\txir`

.. ofile:field:: xi_h
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\Rstar`

   Horizontal displacement perturbation :math:`\txih`

.. ofile:field:: eul_Phi
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar/\Rstar`

   Eulerian potential perturbation :math:`\tPhi'`

.. ofile:field:: deul_Phi
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar/\Rstar^{2}`

   Eulerian potential gradient perturbation :math:`\sderiv{\tPhi'}{x}`

.. ofile:field:: lag_S
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\cP`

   Lagrangian specific entropy perturbation :math:`\delta\tS`

.. ofile:field:: lag_L
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\Lstar`

   Lagrangian radiative luminosity perturbation :math:`\delta\tLrad`

.. ofile:field:: eul_P
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`P`

   Eulerian total pressure perturbation :math:`\tP'`

.. ofile:field:: eul_rho
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\rho`

   Eulerian density perturbation :math:`\trho'`

.. ofile:field:: eul_T
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`T`

   Eulerian temperature perturbation :math:`\tT'`

.. ofile:field:: lag_P
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`P`

   Lagrangian total pressure perturbation :math:`\delta\tP`

.. ofile:field:: lag_rho
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`\rho`

   Lagrangian density perturbation :math:`\delta\trho`

.. ofile:field:: lag_T
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`T`

   Lagrangian temperature perturbation :math:`\delta\tT`

Energetics & Transport
----------------------

.. ofile:field:: eta
   :type: real

   Normalized growth rate :math:`\eta`; evaluated using expression in
   text of page 1186 of :ads_citet:`stellingwerf:1978`\ [#only-N]_

.. ofile:field:: E
   :type: real
   :units: :math:`\Mstar\Rstar^{2}`

   Mode inertia :math:`E \equiv \int \sderiv{E}{x} \, \diff{x}`

.. ofile:field:: E_p
   :type: real
   :units: :math:`\Mstar\Rstar^{2}`

   Acoustic mode inertia :math:`E_{\rm p}`; evaluated by integrating
   :math:`\sderiv{E}{x}` with respect to :math:`x` across regions
   where :math:`\varpi=1`

.. ofile:field:: E_g
   :type: real
   :units: :math:`\Mstar\Rstar^{2}`

   Gravity mode inertia :math:`E_{\rm g}`; evaluated by integrating
   :math:`\sderiv{E}{x}` with respect to :math:`x` across regions
   where :math:`\varpi=-1`

.. ofile:field:: E_norm
   :type: real

   Normalized inertia :math:`E_{\rm norm}`; evaluation controlled by
   :nml:option:`inertia_norm <osc.inertia_norm>` option

.. ofile:field:: E_ratio
   :type: real

   Ratio of mode inertia outside reference location, to total inertia

.. ofile:field:: H
   :type: real
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode energy :math:`H \equiv \frac{1}{2} \omega^{2} E`

.. ofile:field:: W
   :type: real
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode work :math:`W \equiv \int \sderiv{W}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: W_eps
   :type: real
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode nuclear work :math:`W_{\epsilon} \equiv \int
   \sderiv{W_{\epsilon}}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: tau_ss
   :type: real
   :units: :math:`G\Mstar^{2}/\Rstar`

   Steady-state torque :math:`\tau_{\rm ss} \equiv \int
   \sderiv{\tau_{\rm ss}}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: tau_tr
   :type: real
   :units: :math:`G\Mstar^{2}/\Rstar`

   Transient torque :math:`\tau_{\rm tr} \equiv \int
   \sderiv{\tau_{\rm tr}}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: dE_dx
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\Mstar \Rstar^{2}`

   Differential inertia :math:`\sderiv{E}{x}`; evaluated using eqn. (3.139) of
   :ads_citet:`aerts:2010`

.. ofile:field:: dW_dx
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Differential work :math:`\sderiv{W}{x}`; evaluated using
   eqn. (25.9) of :ads_citet:`unno:1989`\ [#only-N]_

.. ofile:field:: dW_eps_dx
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Differential nuclear work :math:`\sderiv{W_{\epsilon}}{x}`;
   evaluated using eqn. (25.9) of :ads_citet:`unno:1989`\ [#only-N]_

.. ofile:field:: dtau_ss_dx
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Steady-state differential torque :math:`\sderiv{\tau_{\rm ss}}{x}`;
   evaluated using eqn. (13) of :ads_citet:`townsend:2018`\ [#only-N]_

.. ofile:field:: dtau_tr_dx
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Transient differential torque :math:`\sderiv{\tau_{\rm tr}}{x}`;
   evaluated using eqn. (14) of :ads_citet:`townsend:2018`\ [#only-N]_

.. ofile:field:: alpha_0
   :type: real
   :dim: :ofile:field:`n`

   Excitation coefficient :math:`\alpha_{0}`; evaluated using
   eqn. (26.10) of :ads_citet:`unno:1989`\ [#only-N]_

.. ofile:field:: alpha_1
   :type: real
   :dim: :ofile:field:`n`

   Excitation coefficient :math:`\alpha_{1}`; evaluated using eqn. (26.12) of
   :ads_citet:`unno:1989`\ [#only-N]_

Rotation
--------

.. ofile:field:: Omega_rot_ref
   :type: real
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Rotation angular frequency :math:`\Omega_{\rm rot,ref}` at
   reference location

.. ofile:field:: Omega_rot
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Rotation angular frequency :math:`\Orot`

.. ofile:field:: domega_rot
   :type: real

   Dimensionless first-order rotational splitting :math:`\Delta
   \omega`; evaluated using eqn. (3.355) of :ads_citet:`aerts:2010`

.. ofile:field:: dfreq_rot
   :type: real
   :units: controlled by :nml:option:`freq_units` and :nml:option:`freq_frame` options

   Dimensioned first-order rotational splitting

.. ofile:field:: beta
   :type: real

   Rotation splitting coefficient :math:`\beta \equiv \int \sderiv{\beta}{x} \, \diff{x}`

.. ofile:field:: dbeta_dx
   :type: real
   :dim: :ofile:field:`n`

   Un-normalized rotation splitting kernel :math:`\sderiv{\beta}{x}`;
   evaluated using eqn. (3.357) of :ads_citet:`aerts:2010`

.. ofile:field:: lambda
   :type: complex
   :dim: :ofile:field:`n`

   Angular eigenvalue :math:`\lambda`

Stellar Structure
-----------------

.. ofile:field:: M_star
   :type: real
   :units: :math:`\gram`

   Stellar mass :math:`\Mstar`\ [#only-D]_

.. ofile:field:: R_star
   :type: real
   :units: :math:`\cm`

   Stellar radius :math:`\Rstar`\ [#only-D]_

.. ofile:field:: L_star
   :type: real
   :units: :math:`\erg\,\second^{-1}`

   Stellar luminosity :math:`\Lstar`\ [#only-D]_

.. ofile:field:: Delta_p
   :type: real
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Asymptotic p-mode large frequency separation :math:`\Delta \nu`

.. ofile:field:: Delta_g
   :type: real
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Asymptotic g-mode inverse period separation :math:`(\Delta P)^{-1}`

.. ofile:field:: V_2
   :type: real
   :dim: :ofile:field:`n`

   Rescaled homology invariant :math:`V_2 \equiv x^{-2} V`; defined in
   :ref:`osc-struct-coeffs` section

.. ofile:field:: As
   :type: real
   :dim: :ofile:field:`n`

   Schwarzschild discriminant :math:`A^{*}`; defined in
   :ref:`osc-struct-coeffs` section

.. ofile:field:: U
   :type: real
   :dim: :ofile:field:`n`

   Homology invariant :math:`U`; defined in :ref:`osc-struct-coeffs`
   section

.. ofile:field:: c_1
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`c_{1}`; defined in
   :ref:`osc-struct-coeffs` section

.. ofile:field:: U_D
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`UD \equiv U \sderiv{\ln\rho}{\ln x}`\ [#only-P]_

.. ofile:field:: Gamma_1
   :type: real
   :dim: :ofile:field:`n`

   Adiabatic exponent :math:`\Gammi`; defined in
   :ref:`osc-linear-eqns` section

.. ofile:field:: upsilon_T
   :type: real
   :dim: :ofile:field:`n`

   Thermodynamic coefficient :math:`\upsT`; defined in
   :ref:`osc-linear-eqns` section\ [#only-N]_

.. ofile:field:: nabla_ad
   :type: real
   :dim: :ofile:field:`n`

   Adiabatic temperature gradient :math:`\nabad`; defined in
   :ref:`osc-linear-eqns` section\ [#only-N]_

.. ofile:field:: dnabla_ad
   :type: real
   :dim: :ofile:field:`n`

   Logarithmic gradient :math:`\dnabad \equiv \sderiv{\ln\nabad}{\ln
   x}` of adiabatic temperature gradient\ [#only-N]_

.. ofile:field:: beta_rad
   :type: real
   :dim: :ofile:field:`n`

   Ratio of radiation pressure to gas pressure [#only-D]_

.. ofile:field:: nabla
   :type: real
   :dim: :ofile:field:`n`

   Temperature gradient :math:`\nabla`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_lum
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\clum`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_rad
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\crad`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_thn
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\cthn`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_thk
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\cthk`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_eps
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\ceps`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: c_egv
   :type: real
   :dim: :ofile:field:`n`

   Structure coefficient :math:`\cegv`; defined in
   :ref:`osc-struct-coeffs` section\ [#only-N]_

.. ofile:field:: eps_rho
   :type: real
   :dim: :ofile:field:`n`

   Nuclear energy generation partial :math:`\epsnucrho`; defined in
   :ref:`osc-linear-eqns` section\ [#only-N]_

.. ofile:field:: eps_T
   :type: real
   :dim: :ofile:field:`n`

   Nuclear energy generation partial :math:`\epsnucT`; defined in
   :ref:`osc-linear-eqns` section\ [#only-N]_

.. ofile:field:: kap_rho
   :type: real
   :dim: :ofile:field:`n`

   Opacity partial :math:`\kaprho`; defined in :ref:`osc-linear-eqns`
   section\ [#only-N]_

.. ofile:field:: kap_T
   :type: real
   :dim: :ofile:field:`n`

   Opacity partial :math:`\kapT`; defined in :ref:`osc-linear-eqns`
   section\ [#only-N]_

.. ofile:field:: M_r
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\gram`

   Interior mass :math:`M_r`\ [#only-D]_

.. ofile:field:: P
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\barye`

   Total pressure :math:`P`\ [#only-D]_

.. ofile:field:: rho
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\gram\,\cm^{-3}`

   Density :math:`\rho`\ [#only-D]_

.. ofile:field:: T
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`\kelvin`

   Temperature :math:`T`\ [#only-D]_

Tidal Response
--------------

Note that these fields are available only when using :program:`gyre_tides`.

.. ofile:field:: k
   :type: integer

   Fourier harmonic :math:`k`

.. ofile:field:: eul_Psi_ref
   :type: complex
   :units: :math:`G\Mstar/\Rstar`

   Eulerian total potential perturbation :math:`\tPsi'_{\rm ref}` at
   reference location

.. ofile:field:: Phi_T_ref
   :type: real
   :units: :math:`G\Mstar/\Rstar`

   Tidal potential :math:`\tPhi_{\rm T, ref}` at reference location

.. ofile:field:: eul_Psi
   :type: complex
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar/\Rstar`

   Eulerian total potential perturbation :math:`\tPsi'`

.. ofile:field:: Phi_T
   :type: real
   :dim: :ofile:field:`n`
   :units: :math:`G\Mstar/\Rstar`

   Tidal potential :math:`\tPhi_{{\rm T}}`

.. ofile:field:: Omega_orb
   :type: real
   :units: controlled by :nml:option:`freq_units` and :nml:option:`freq_frame` options

   Orbital angular frequency :math:`\Oorb`

.. ofile:field:: q
   :type: real

   Ratio :math:`q` of secondary mass to primary mass

.. ofile:field:: e
   :type: real

   Orbital eccentricity :math:`e`

.. ofile:field:: R_a
   :type: real

   Ratio :math:`R/a` of primary radius to orbital semi-major axis

.. ofile:field:: cbar
   :type: real

   Tidal expansion coefficient :math:`\cbar_{\ell,m,k}`; see eqn. (A1) of :ads_citet:`sun:2023`

.. ofile:field:: Gbar_1
   :type: real

   Secular orbital evolution coefficient
   :math:`\Gbar^{(1)}_{\ell,m,k}`; equivalent to
   :math:`G^{(1)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_2
   :type: real

   Secular orbital evolution coefficient
   :math:`\Gbar^{(2)}_{\ell,m,k}`; equivalent to
   :math:`G^{(2)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_3
   :type: real

   Secular orbital evolution coefficient
   :math:`\Gbar^{(3)}_{\ell,m,k}`; equivalent to
   :math:`G^{(3)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_4
   :type: real

   Secular orbital evolution coefficient
   :math:`\Gbar^{(4)}_{\ell,m,k}`; equivalent to
   :math:`G^{(4)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_E
   :type: real

   Secular energy transfer coefficient :math:`\Gbar^{(E)}_{\ell,m,k}`;
   derived from :math:`\Gbar^{(4)}_{\ell,m,k}` by dropping the leading
   :math:`m` term

.. rubric:: Footnotes

.. [#only-N] This field is available only for stellar models with
             :ref:`N capability <model-caps>`

.. [#only-D] This field is available only for stellar models with
             :ref:`D capability <model-caps>`

.. [#only-P] This field is available only for :ref:`homogeneous
             <hom-models>`, :ref:`polytropic <poly-models>` and
             :ref:`analytic polytropic <anapoly-models>` models
