.. _summary-files:

.. ofile:filetype:: summary

.. nml:group:: output
   :no-target:

Summary Files
=============

*Summary* files collect together global properties, such as
eigenfrequencies and radial orders, of all solutions (modes, tidal
responses, etc.) found during a run. The specific data written to a
summary file are controlled by the :nml:option:`summary_item_list`
options of the :nml:group:`ad_output` and :nml:group:`nad_output` namelist
groups (:program:`gyre` adiabatic and non-adiabatic calculations,
respectively) and the :nml:group:`tide_output` namelist group
(:program:`gyre_tides` calculations). These options specify the
data to be written via a comma-separated list of fields.

The following subsections describe the fields that may appear in
:nml:option:`summary_item_list`, grouped together by functional area.

Solution Data
-------------

.. ofile:field:: n_row
   :type: integer

   Number of rows :math:`N_{\rm row}` in summary file, each
   corresponding to a mode found (:program:`gyre`) or a tidal
   response evaluated (:program:`gyre_tides`)

.. ofile:field:: n
   :type: integer
   :dim: :ofile:field:`n_row`

   Number of spatial grid points :math:`N`

.. ofile:field:: omega
   :type: complex
   :dim: :ofile:field:`n_row`

   Dimensionless angular frequency :math:`\omega`

.. ofile:field:: x_ref
   :type: real
   :dim: :ofile:field:`n_row`

   Dimensionless radial coordinate :math:`\xref` of reference location

.. ofile:field:: chi
   :type: real
   :dim: :ofile:field:`n_row`

   Root-finding convergence parameter :math:`\chi`

.. ofile:field:: n_iter
   :type: integer
   :dim: :ofile:field:`n_row`

   Root-finding number of iterations

Observables
-----------

.. ofile:field:: freq
   :type: complex
   :dim: :ofile:field:`n_row`
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
   :dim: :ofile:field:`n_row`

   Effective temperature perturbation amplitude :math:`f_{T}`;
   evaluated at reference location using eqn. (5) of
   :ads_citet:`dupret:2003`

.. ofile:field:: f_g
   :type: real
   :dim: :ofile:field:`n_row`

   Effective gravity perturbation amplitude :math:`f_{g}`; evaluated
   at reference location using eqn. (6) of :ads_citet:`dupret:2003`

.. ofile:field:: psi_T
   :type: real
   :dim: :ofile:field:`n_row`

   Effective temperature perturbation phase :math:`\psi_{T}`;
   evaluated at reference location using eqn. (5) of
   :ads_citet:`dupret:2003`

.. ofile:field:: psi_g
   :type: real
   :dim: :ofile:field:`n_row`

   Effective gravity perturbation phase :math:`\psi_{g}`; evaluated at
   reference location using eqn. (6) of :ads_citet:`dupret:2003`

Classification & Validation
---------------------------

.. ofile:field:: id
   :type: integer
   :dim: :ofile:field:`n_row`

   Unique mode index

.. ofile:field:: l
   :type: integer
   :dim: :ofile:field:`n_row`

   Harmonic degree :math:`\ell`

.. ofile:field:: l_i
   :type: complex
   :dim: :ofile:field:`n_row`

   Effective harmonic degree at inner boundary :math:`\ell_{\rm i}`

.. ofile:field:: m
   :type: integer
   :dim: :ofile:field:`n_row`

   Azimuthal order :math:`m`

.. ofile:field:: n_p
   :type: integer
   :dim: :ofile:field:`n_row`

   Acoustic-wave winding number :math:`\nump`

.. ofile:field:: n_g
   :type: integer
   :dim: :ofile:field:`n_row`

   Gravity-wave winding number :math:`\numg`

.. ofile:field:: n_pg
   :type: integer
   :dim: :ofile:field:`n_row`

   Radial order :math:`\numpg` within the
   Eckart-Scuflaire-Osaki-Takata scheme (see
   :ads_citealp:`takata:2006b`)

.. ofile:field:: omega_int
   :type: complex
   :dim: :ofile:field:`n_row`

   Dimensionless eigenfrequency :math:`\omega_{\rm int}` based on
   integral expression; evaluated using eqn. (A8) of
   :ads_citet:`townsend:2025`

.. ofile:field:: zeta
   :type: complex
   :dim: :ofile:field:`n_row`

   Integrated frequency weight :math:`\zeta \equiv \int \sderiv{\zeta}{x} \, \diff{x}`

Perturbations
-------------

.. ofile:field:: xi_r_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`\Rstar`

   Radial displacement perturbation :math:`\txi_{r,{\rm ref}}` at
   reference location

.. ofile:field:: xi_h_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`\Rstar`

   Horizontal displacement perturbation :math:`\txi_{\rm h,ref}` at
   reference location

.. ofile:field:: eul_Phi_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units:  :math:`G\Mstar/\Rstar`

   Eulerian potential perturbation :math:`\tPhi'_{\rm ref}` at
   reference location

.. ofile:field:: deul_Phi_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar/\Rstar^{2}`

   Eulerian potential gradient perturbation
   :math:`(\sderiv{\tPhi'}{x})_{\rm ref}` at reference location

.. ofile:field:: lag_S_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`\cP`

   Lagrangian specific entropy perturbation :math:`\delta\tS_{\rm
   ref}` at reference location

.. ofile:field:: lag_L_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`\Lstar`

   Lagrangian radiative luminosity perturbation :math:`\delta\tL_{\rm
   R,ref}` at reference location

Energetics & Transport
----------------------

.. ofile:field:: eta
   :type: real
   :dim: :ofile:field:`n_row`

   Normalized growth rate :math:`\eta`; evaluated using expression in
   text of page 1186 of :ads_citet:`stellingwerf:1978`\ [#only-N]_

.. ofile:field:: E
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\Mstar\Rstar^{2}`

   Mode inertia :math:`E \equiv \int \sderiv{E}{x} \, \diff{x}`

.. ofile:field:: E_p
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\Mstar\Rstar^{2}`

   Acoustic mode inertia :math:`E_{\rm p}`; evaluated by integrating
   :math:`\sderiv{E}{x}` with respect to :math:`x` across regions
   where :math:`\varpi=1`

.. ofile:field:: E_g
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\Mstar\Rstar^{2}`

   Gravity mode inertia :math:`E_{\rm g}`; evaluated by integrating
   :math:`\sderiv{E}{x}` with respect to :math:`x` across regions
   where :math:`\varpi=-1`

.. ofile:field:: E_norm
   :type: real
   :dim: :ofile:field:`n_row`

   Normalized inertia :math:`E_{\rm norm}`; evaluation controlled by
   :nml:option:`inertia_norm <osc.inertia_norm>` option

.. ofile:field:: E_ratio
   :type: real
   :dim: :ofile:field:`n_row`

   Ratio of mode inertia outside reference location, to total inertia

.. ofile:field:: H
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode energy :math:`H \equiv \frac{1}{2} \omega^{2} E`

.. ofile:field:: W
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode work :math:`W \equiv \int \sderiv{W}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: W_eps
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Mode nuclear work :math:`W_{\epsilon} \equiv \int
   \sderiv{W_{\epsilon}}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: tau_ss
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar^{2}/\Rstar`

   Steady-state torque :math:`\tau_{\rm ss} \equiv \int
   \sderiv{\tau_{\rm ss}}{x} \, \diff{x}`\ [#only-N]_

.. ofile:field:: tau_tr
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar^{2}/\Rstar`

Rotation
--------

.. ofile:field:: Omega_rot_ref
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Rotation angular frequency :math:`\Omega_{\rm rot,ref}` at
   reference location

.. ofile:field:: domega_rot
   :type: real
   :dim: :ofile:field:`n_row`

   Dimensionless first-order rotational splitting :math:`\Delta
   \omega`; evaluated using eqn. (3.355) of :ads_citet:`aerts:2010`

.. ofile:field:: dfreq_rot
   :type: real
   :dim: :ofile:field:`n_row`
   :units: controlled by :nml:option:`freq_units` and :nml:option:`freq_frame` options

   Dimensioned first-order rotational splitting

.. ofile:field:: beta
   :type: real
   :dim: :ofile:field:`n_row`

   Rotation splitting coefficient :math:`\beta \equiv \int \sderiv{\beta}{x} \, \diff{x}`

Stellar Structure
-----------------

.. ofile:field:: M_star
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\gram`

   Stellar mass :math:`\Mstar`\ [#only-D]_

.. ofile:field:: R_star
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\cm`

   Stellar radius :math:`\Rstar`\ [#only-D]_

.. ofile:field:: L_star
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\erg\,\second^{-1}`

   Stellar luminosity :math:`\Lstar`\ [#only-D]_

.. ofile:field:: Delta_p
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Asymptotic p-mode large frequency separation :math:`\Delta \nu`

.. ofile:field:: Delta_g
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`\sqrt{G\Mstar/\Rstar^{3}}`

   Asymptotic g-mode inverse period separation :math:`(\Delta P)^{-1}`

Tidal Response
--------------

Note that these fields are available only when using :program:`gyre_tides`.

.. ofile:field:: k
   :type: integer
   :dim: :ofile:field:`n_row`

   Fourier harmonic :math:`k`

.. ofile:field:: eul_Psi_ref
   :type: complex
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar/\Rstar`

   Eulerian total potential perturbation :math:`\tPsi'_{\rm ref}` at
   reference location

.. ofile:field:: Phi_T_ref
   :type: real
   :dim: :ofile:field:`n_row`
   :units: :math:`G\Mstar/\Rstar`

   Tidal potential :math:`\tPhi_{\rm T, ref}` at reference location

.. ofile:field:: Omega_orb
   :type: real
   :dim: :ofile:field:`n_row`
   :units: controlled by :nml:option:`freq_units` and :nml:option:`freq_frame` options

   Orbital angular frequency :math:`\Oorb`

.. ofile:field:: q
   :type: real
   :dim: :ofile:field:`n_row`

   Ratio :math:`q` of secondary mass to primary mass

.. ofile:field:: e
   :type: real
   :dim: :ofile:field:`n_row`

   Orbital eccentricity :math:`e`

.. ofile:field:: R_a
   :type: real
   :dim: :ofile:field:`n_row`

   Ratio :math:`R/a` of primary radius to orbital semi-major axis

.. ofile:field:: cbar
   :type: real
   :dim: :ofile:field:`n_row`

   Tidal expansion coefficient :math:`\cbar_{\ell,m,k}`; see eqn. (A1) of :ads_citet:`sun:2023`

.. ofile:field:: Gbar_1
   :type: real
   :dim: :ofile:field:`n_row`

   Secular orbital evolution coefficient
   :math:`\Gbar^{(1)}_{\ell,m,k}`; equivalent to
   :math:`G^{(1)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_2
   :type: real
   :dim: :ofile:field:`n_row`

   Secular orbital evolution coefficient
   :math:`\Gbar^{(2)}_{\ell,m,k}`; equivalent to
   :math:`G^{(2)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_3
   :type: real
   :dim: :ofile:field:`n_row`

   Secular orbital evolution coefficient
   :math:`\Gbar^{(3)}_{\ell,m,k}`; equivalent to
   :math:`G^{(3)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_4
   :type: real
   :dim: :ofile:field:`n_row`

   Secular orbital evolution coefficient
   :math:`\Gbar^{(4)}_{\ell,m,k}`; equivalent to
   :math:`G^{(4)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. ofile:field:: Gbar_E
   :type: real
   :dim: :ofile:field:`n_row`

   Secular energy transfer coefficient :math:`\Gbar^{(E)}_{\ell,m,k}`;
   derived from :math:`\Gbar^{(4)}_{\ell,m,k}` by dropping the leading
   :math:`m` term

.. rubric:: Footnotes

.. [#only-N] This field is available only for stellar models with
             :ref:`N capability <model-caps>`

.. [#only-D] This field is available only for stellar models with
             :ref:`D capability <model-caps>`
