.. _summary-files:

Summary Files
=============

The data written to summary files are controlled by the
:nml_n:`summary_item_list` parameter of the :nml_g:`ad_output`
namelist group (for adiabatic calculations) and the
:nml_g:`nad_output` namelist group (for nonadiabatic
calculations). This parameter is a comma-separated list of items to
appear in the summary file; the following subsections describe the
items that may appear, grouped together by functional area. For each
item, the corresponding math symbol is given (if there is one),
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
   * - :nml_v:`n_j`
     - :math:`N_{j}`
     - integer
     - number of modes found
   * - :nml_v:`omega`
     - :math:`\omega`
     - complex(:nml_v:`n_j`)
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
     - complex(:nml_v:`n_j`)
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
     - real(:nml_v:`n_j`)
     - Effective temperature perturbation amplitude; evaluated using
       eqn. 5 of :ads_citet:`dupret:2003`
   * - :nml_v:`f_g`
     - :math:`f_{\rm g}`
     - real(:nml_v:`n_j`)
     - Effective gravity perturbation amplitude; evaluated using
       eqn. 6 of :ads_citet:`dupret:2003`
   * - :nml_v:`\psi_T`
     - :math:`\psi_{T}`
     - real(:nml_v:`n_j`)
     - Effective temperature perturbation phase; evaluated using
       eqn. 5 of :ads_citet:`dupret:2003`
   * - :nml_v:`f_g`
     - :math:`\psi_{\rm g}`
     - real(:nml_v:`n_j`)
     - Effective gravity perturbation phase; evaluated using
       eqn. 6 of :ads_citet:`dupret:2003`

Classification & Validation
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 65

   * - :nml_v:`j`
     - :math:`j`
     - integer(:nml_v:`n_j`)
     - unique mode index
   * - :nml_v:`l`
     - :math:`\ell`
     - integer(:nml_v:`n_j`)
     - harmonic degree
   * - :nml_v:`l_i`
     - :math:`\ell_{\rm i}`
     - complex(:nml_v:`n_j`)
     - effective harmonic degree at inner boundary
   * - :nml_v:`m`
     - :math:`m`
     - integer(:nml_v:`n_j`)
     - azimuthal order
   * - :nml_v:`n_p`
     - :math:`\np`
     - integer(:nml_v:`n_j`)
     - acoustic-wave winding number
   * - :nml_v:`n_g`
     - :math:`\ng`
     - integer(:nml_v:`n_j`)
     - gravity-wave winding number
   * - :nml_v:`n_pg`
     - :math:`\npg`
     - integer(:nml_v:`n_j`)
     - radial order within the Eckart-Scuflaire-Osaki-Takata
       scheme (see :ads_citealp:`takata:2006b`)
   * - :nml_v:`omega_int`
     - :math:`\omega_{\rm int}`
     - complex(:nml_v:`n_j`)
     - dimensionless eigenfrequency; evaluated by
       integrating :math:`\sderiv{\zeta}{x}`

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
     - complex(:nml_v:`n_j`)
     - radial displacement perturbation at reference location [:math:`R`]
   * - :nml_v:`eul_phi_ref`
     - :math:`\tPhi'_{\rm ref}`
     - complex(:nml_v:`n_j`)
     - Eulerian potential perturbation at reference location [:math:`GM/R`]
   * - :nml_v:`deul_phi_ref`
     - :math:`(\sderiv{\tPhi'}{x})_{\rm ref}`
     - complex(:nml_v:`n_j`)
     - Eulerian potential gradient perturbation at reference location [:math:`GM/R^{2}`]
   * - :nml_v:`lag_S_ref`
     - :math:`\delta\tS_{\rm ref}`
     - complex(:nml_v:`n_j`)
     - Lagrangian specific entropy perturbation at reference location [:math:`R`]
   * - :nml_v:`lag_L_ref`
     - :math:`\delta\tL_{\rm R,ref}`
     - complex(:nml_v:`n_j`)
     - Lagrangian radiative luminosity perturbation at reference location [:math:`L`]

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
     - real(:nml_v:`n_j`)
     - normalized growth rate :math:`\eta`; evaluated using expression
       in text of page 1186 of :ads_citet:`stellingwerf:1978`
   * - :nml_v:`E`
     - :math:`E`
     - real(:nml_v:`n_j`)
     - mode inertia [:math:`M R^{2}`]; evaluated by integrating
       :math:`\sderiv{E}{x}`
   * - :nml_v:`E_p`
     - :math:`E_{\rm p}`
     - real(:nml_v:`n_j`)
     - acoustic mode inertia [:math:`M R^{2}`]; evaluated by
       integrating :math:`\sderiv{E}{x}` where
       :math:`\varpi=1`
   * - :nml_v:`E_g`
     - :math:`E_{\rm g}`
     - real(:nml_v:`n_j`)
     - gravity mode inertia [:math:`M R^{2}`]; evaluated by
       integrating :math:`\sderiv{E}{x}` in regions wherre
       :math:`\varpi=-1`
   * - :nml_v:`E_norm`
     - :math:`E_{\rm norm}`
     - real(:nml_v:`n_j`)
     - normalized inertia; evaluation controlled by :nml_n:`inertia_norm`
       parameter
   * - :nml_v:`E_ratio`
     - ---
     - real(:nml_v:`n_j`)
     - ratio of mode inertias inertia inside/outside reference
       location
   * - :nml_v:`H`
     - :math:`H`
     - real(:nml_v:`n_j`)
     - mode energy [:math:`G M^{2}/R`]; evaluated as
       :math:`\frac{1}{2} \omega^{2} E`
   * - :nml_v:`W`\ [#only-N]_
     - :math:`W`
     - real(:nml_v:`n_j`)
     - mode work [:math:`G M^{2}/R`]; evaluated by
       integrating :math:`\sderiv{W}{x}`
   * - :nml_v:`W_eps`\ [#only-N]_
     - :math:`W_{\epsilon}`
     - real(:nml_v:`n_j`)
     - mode work [:math:`G M^{2}/R`]; evaluated by
       integrating :math:`\sderiv{W_{\epsilon}}{x}`
   * - :nml_v:`tau_ss`
     - :math:`\tau_{\rm ss}`
     - real(:nml_v:`n_j`)
     - steady-state torque [:math:`G M^{2}/R`]; evaluated by
       integrating :math:`\sderiv{\tau_{\rm ss}}{x}`
   * - :nml_v:`tau_tr`
     - :math:`\tau_{\rm tr}`
     - real(:nml_v:`n_j`)
     - steady-state torque [:math:`G M^{2}/R`]; evaluated by
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
   * - :nml_v:`domega_rot`
     - :math:`\delta \omega`
     - real(:nml_v:`n_j`)
     - dimensionless first-order rotational splitting; evaluated using eqn. 3.355 of :ads_citet:`aerts:2010`
   * - :nml_v:`dfreq_rot`
     - ---
     - real(:nml_v:`n_j`)
     - dimensioned first-order rotational splitting; units and reference frame controlled by
       :nml_n:`freq_units` and :nml_n:`freq_frame` parameters
   * - :nml_v:`beta`
     - :math:`\beta`
     - real(:nml_v:`n_j`)
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
     - :math:`M`
     - real(:nml_v:`n_j`)
     - stellar mass [:math:`\gram`]
   * - :nml_v:`R_star`\ [#only-D]_
     - :math:`R`
     - real(:nml_v:`n_j`)
     - stellar radiua [:math:`\cm`]
   * - :nml_v:`L_star`\ [#only-D]_
     - :math:`L`
     - real(:nml_v:`n_j`)
     - stellar luminosity [:math:`\erg\,\second^{-1}`]
   * - :nml_v:`Delta_p`
     - :math:`\Delta \nu`
     - real(:nml_v:`n_j`)
     - asymptotic p-mode large frequency separation [:math:`\sqrt{GM/R^{3}}`]
   * - :nml_v:`Delta_g`
     - :math:`(\Delta P)^{-1}`
     - real(:nml_v:`n_j`)
     - asymptotic g-mode inverse period separation [:math:`\sqrt{GM/R^{3}}`]

.. rubric:: Footnotes

.. [#only-N] This option is available only for stellar models with :ref:`N capability <model-caps>`

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
