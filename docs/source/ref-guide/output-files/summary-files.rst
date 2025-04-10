.. _summary-files:

Summary Files
=============

*Summary* files collect together global properties, such as
eigenfrequencies and radial orders, of all solutions (modes, tidal
responses, etc.) found during a run. The specific data written to a
summary file are controlled by the :nml_n:`summary_item_list`
parameters of the :nml_g:`ad_output` and :nml_g:`nad_output` namelist
groups (:program:`gyre` adiabatic and non-adiabatic calculations,
respectively) and the :nml_g:`tides_output` namelist group
(:program:`gyre_tides` calculations). These parameters specify the
items to be written, via a comma-separated list.

The following subsections describe the items that may appear in
:nml_n:`summary_item_list`, grouped together by functional area. For
each item, the corresponding math symbol is given (if there is one),
together with the datatype, and a brief description. Units (where
applicable) are indicated in brackets [].

Solution Data
-------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`n_row`
     - :math:`N_{\rm row}`
     - integer
     - number of rows in summary file, each corresponding to a mode found
       (:program:`gyre`) or a tidal response evaluated (:program:`gyre_tides`)
   * - :nml_v:`n`
     - :math:`N`
     - integer(:nml_v:`n_row`)
     - number of spatial grid points
   * - :nml_v:`omega`
     - :math:`\omega`
     - complex(:nml_v:`n_row`)
     - dimensionless eigenfrequency

Observables
-----------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`freq`
     - ---
     - complex(:nml_v:`n_row`)
     - dimensioned frequency; units and reference frame controlled by
       :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
   * - :nml_v:`freq_units`
     - ---
     - string
     - :nml_n:`freq_units` parameter
   * - :nml_v:`freq_frame`
     - ---
     - string
     - :nml_n:`freq_frame` parameter
   * - :nml_v:`f_T`
     - :math:`f_{T}`
     - real(:nml_v:`n_row`)
     - Effective temperature perturbation amplitude; evaluated using
       eqn. 5 of :ads_citet:`dupret:2003`
   * - :nml_v:`f_g`
     - :math:`f_{\rm g}`
     - real(:nml_v:`n_row`)
     - Effective gravity perturbation amplitude; evaluated using
       eqn. 6 of :ads_citet:`dupret:2003`
   * - :nml_v:`psi_T`
     - :math:`\psi_{T}`
     - real(:nml_v:`n_row`)
     - Effective temperature perturbation phase; evaluated using
       eqn. 5 of :ads_citet:`dupret:2003`
   * - :nml_v:`psi_g`
     - :math:`\psi_{\rm g}`
     - real(:nml_v:`n_row`)
     - Effective gravity perturbation phase; evaluated using
       eqn. 6 of :ads_citet:`dupret:2003`

Classification & Validation
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`id`
     - ---
     - integer(:nml_v:`n_row`)
     - unique mode index
   * - :nml_v:`l`
     - :math:`\ell`
     - integer(:nml_v:`n_row`)
     - harmonic degree
   * - :nml_v:`l_i`
     - :math:`\ell_{\rm i}`
     - complex(:nml_v:`n_row`)
     - effective harmonic degree at inner boundary
   * - :nml_v:`m`
     - :math:`m`
     - integer(:nml_v:`n_row`)
     - azimuthal order
   * - :nml_v:`n_p`
     - :math:`\nump`
     - integer(:nml_v:`n_row`)
     - acoustic-wave winding number
   * - :nml_v:`n_g`
     - :math:`\numg`
     - integer(:nml_v:`n_row`)
     - gravity-wave winding number
   * - :nml_v:`n_pg`
     - :math:`\numpg`
     - integer(:nml_v:`n_row`)
     - radial order within the Eckart-Scuflaire-Osaki-Takata
       scheme (see :ads_citealp:`takata:2006b`)
   * - :nml_v:`omega_int`
     - :math:`\omega_{\rm int}`
     - complex(:nml_v:`n_row`)
     - dimensionless eigenfrequency based on integral expression; evaluated using eqn. A8 of Townsend et al. (2025)
   * - :nml_v:`zeta`
     - :math:`\zeta`
     - complex(:nml_v:`n_row`)
     - integral of :math:`\sderiv{\zeta}{x}` with respect to :math:`x`

Perturbations
-------------
  
.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`x_ref`
     - :math:`x_{\rm ref}`
     - real
     - fractional radius of reference location
   * - :nml_v:`xi_r_ref`
     - :math:`\txi_{r,{\rm ref}}`
     - complex(:nml_v:`n_row`)
     - radial displacement perturbation at reference location [:math:`\Rstar`]
   * - :nml_v:`xi_h_ref`
     - :math:`\txi_{h,{\rm ref}}`
     - complex(:nml_v:`n_row`)
     - radial displacement perturbation at reference location [:math:`\Rstar`]
   * - :nml_v:`eul_Phi_ref`
     - :math:`\tPhi'_{\rm ref}`
     - complex(:nml_v:`n_row`)
     - Eulerian potential perturbation at reference location [:math:`G\Mstar/\Rstar`]
   * - :nml_v:`deul_Phi_ref`
     - :math:`(\sderiv{\tPhi'}{x})_{\rm ref}`
     - complex(:nml_v:`n_row`)
     - Eulerian potential gradient perturbation at reference location [:math:`G\Mstar/\Rstar^{2}`]
   * - :nml_v:`lag_S_ref`
     - :math:`\delta\tS_{\rm ref}`
     - complex(:nml_v:`n_row`)
     - Lagrangian specific entropy perturbation at reference location [:math:`\cP`]
   * - :nml_v:`lag_L_ref`
     - :math:`\delta\tL_{\rm R,ref}`
     - complex(:nml_v:`n_row`)
     - Lagrangian radiative luminosity perturbation at reference location [:math:`\Lstar`]

Energetics & Transport
----------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`eta`\ [#only-N]_
     - :math:`\eta`
     - real(:nml_v:`n_row`)
     - normalized growth rate :math:`\eta`; evaluated using expression
       in text of page 1186 of :ads_citet:`stellingwerf:1978`
   * - :nml_v:`E`
     - :math:`E`
     - real(:nml_v:`n_row`)
     - mode inertia [:math:`\Mstar\Rstar^{2}`]; evaluated by integrating
       :math:`\sderiv{E}{x}`
   * - :nml_v:`E_p`
     - :math:`E_{\rm p}`
     - real(:nml_v:`n_row`)
     - acoustic mode inertia [:math:`\Mstar\Rstar^{2}`]; evaluated by
       integrating :math:`\sderiv{E}{x}` where
       :math:`\varpi=1`
   * - :nml_v:`E_g`
     - :math:`E_{\rm g}`
     - real(:nml_v:`n_row`)
     - gravity mode inertia [:math:`\Mstar\Rstar^{2}`]; evaluated by
       integrating :math:`\sderiv{E}{x}` in regions where
       :math:`\varpi=-1`
   * - :nml_v:`E_norm`
     - :math:`E_{\rm norm}`
     - real(:nml_v:`n_row`)
     - normalized inertia; evaluation controlled by :nml_n:`inertia_norm`
       parameter
   * - :nml_v:`E_ratio`
     - ---
     - real(:nml_v:`n_row`)
     - ratio of mode inertia outside reference location, to total inertia
   * - :nml_v:`H`
     - :math:`H`
     - real(:nml_v:`n_row`)
     - mode energy [:math:`G\Mstar^{2}/\Rstar`]; evaluated as
       :math:`\frac{1}{2} \omega^{2} E`
   * - :nml_v:`W`\ [#only-N]_
     - :math:`W`
     - real(:nml_v:`n_row`)
     - mode work [:math:`G\Mstar^{2}/\Rstar`]; evaluated by
       integrating :math:`\sderiv{W}{x}`
   * - :nml_v:`W_eps`\ [#only-N]_
     - :math:`W_{\epsilon}`
     - real(:nml_v:`n_row`)
     - mode work [:math:`G\Mstar^{2}/\Rstar`]; evaluated by
       integrating :math:`\sderiv{W_{\epsilon}}{x}`
   * - :nml_v:`tau_ss`
     - :math:`\tau_{\rm ss}`
     - real(:nml_v:`n_row`)
     - steady-state torque [:math:`G\Mstar^{2}/\Rstar`]; evaluated by
       integrating :math:`\sderiv{\tau_{\rm ss}}{x}`
   * - :nml_v:`tau_tr`
     - :math:`\tau_{\rm tr}`
     - real(:nml_v:`n_row`)
     - steady-state torque [:math:`G\Mstar^{2}/\Rstar`]; evaluated by
       integrating :math:`\sderiv{\tau_{\rm tr}}{x}`

Rotation
--------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`Omega_rot_ref`
     - :math:`\Omega_{\rm rot,ref}`
     - real(:nml_v:`n_row`)
     - rotation angular frequency at reference location[:math:`\sqrt{G\Mstar/\Rstar^{3}}`]
   * - :nml_v:`domega_rot`
     - :math:`\Delta \omega`
     - real(:nml_v:`n_row`)
     - dimensionless first-order rotational splitting; evaluated using eqn. 3.355 of :ads_citet:`aerts:2010`
   * - :nml_v:`dfreq_rot`
     - ---
     - real(:nml_v:`n_row`)
     - dimensioned first-order rotational splitting; units and reference frame controlled by
       :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
   * - :nml_v:`beta`
     - :math:`\beta`
     - real(:nml_v:`n_row`)
     - rotation splitting coefficient; evaluated by
       integrating :math:`\sderiv{\beta}{x}`

Stellar Structure
-----------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`M_star`\ [#only-D]_
     - :math:`\Mstar`
     - real(:nml_v:`n_row`)
     - stellar mass [:math:`\gram`]
   * - :nml_v:`R_star`\ [#only-D]_
     - :math:`\Rstar`
     - real(:nml_v:`n_row`)
     - stellar radius [:math:`\cm`]
   * - :nml_v:`L_star`\ [#only-D]_
     - :math:`\Lstar`
     - real(:nml_v:`n_row`)
     - stellar luminosity [:math:`\erg\,\second^{-1}`]
   * - :nml_v:`Delta_p`
     - :math:`\Delta \nu`
     - real(:nml_v:`n_row`)
     - asymptotic p-mode large frequency separation [:math:`\sqrt{G\Mstar/\Rstar^{3}}`]
   * - :nml_v:`Delta_g`
     - :math:`(\Delta P)^{-1}`
     - real(:nml_v:`n_row`)
     - asymptotic g-mode inverse period separation [:math:`\sqrt{G\Mstar/\Rstar^{3}}`]

Tidal Response
--------------

Note that these items are available only when using :program:`gyre_tides`.

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - Item
     - Symbol
     - Datatype
     - Description
   * - :nml_v:`k`
     - :math:`k`
     - integer(:nml_v:`n_row`)
     - Fourier harmonic
   * - :nml_v:`eul_Psi_ref`
     - :math:`\tPsi'_{\rm ref}`
     - complex(:nml_v:`n_row`)
     - Eulerian total potential perturbation at reference location [:math:`G\Mstar/\Rstar`]
   * - :nml_v:`Phi_T_ref`
     - :math:`\tPhi_{\rm T, ref}`
     - real(:nml_v:`n_row`)
     - tidal potential at reference location [:math:`G\Mstar/\Rstar`]
   * - :nml_v:`Omega_orb`
     - :math:`\Oorb`
     - real(:nml_v:`n_row`)
     - orbital angular frequency; units and reference frame controlled by
       :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
   * - :nml_v:`q`
     - :math:`q`
     - real(:nml_v:`n_row`)
     - ratio of secondary mass to primary mass
   * - :nml_v:`e`
     - :math:`e`
     - real(:nml_v:`n_row`)
     - orbital eccentricity
   * - :nml_v:`R_a`
     - :math:`R/a`
     - real(:nml_v:`n_row`)
     - ratio of primary radius to orbital semi-major axis 
   * - :nml_v:`cbar`
     - :math:`\cbar_{\ell,m,k}`
     - real(:nml_v:`n_row`)
     - tidal expansion coefficient; see eqn. A1 of :ads_citet:`sun:2023`
   * - :nml_v:`Gbar_1`
     - :math:`\Gbar^{(1)}_{\ell,m,k}`
     - real(:nml_v:`n_row`)
     - secular orbital evolution coefficient; equivalent to :math:`G^{(1)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)
   * - :nml_v:`Gbar_2`
     - :math:`\Gbar^{(2)}_{\ell,m,k}`
     - real(:nml_v:`n_row`)
     - secular orbital evolution coefficient; equivalent to :math:`G^{(2)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)
   * - :nml_v:`Gbar_3`
     - :math:`\Gbar^{(3)}_{\ell,m,k}`
     - real(:nml_v:`n_row`)
     - secular orbital evolution coefficient; equivalent to :math:`G^{(3)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)
   * - :nml_v:`Gbar_4`
     - :math:`\Gbar^{(4)}_{\ell,m,k}`
     - real(:nml_v:`n_row`)
     - secular orbital evolution coefficient; equivalent to :math:`G^{(4)}_{\ell,m,-k}` (see :ads_citealp:`willems:2003`)

.. rubric:: Footnotes

.. [#only-N] This item is available only for stellar models with :ref:`N capability <model-caps>`

.. [#only-D] This item is available only for stellar models with :ref:`D capability <model-caps>`
