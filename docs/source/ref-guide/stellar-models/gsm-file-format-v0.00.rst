Version 0.00
------------

Data items in the root HDF5 group of version-0.00 GSM-format files are as follows:

.. list-table::
   :widths: 10 10 10 10 60
   :header-rows: 1

   * - Item
     - Symbol
     - Object type
     - Data type
     - Definition
   * - :code:`n`
     - :math:`N`
     - attribute
     - integer
     - number of grid points
   * - :code:`M_star`
     - :math:`\Mstar`
     - attribute
     - real
     - stellar mass [:math:`\gram`]
   * - :code:`R_star`
     - :math:`\Rstar`
     - attribute
     - real
     - photospheric radius [:math:`\cm`]
   * - :code:`L_star`
     - :math:`\Lstar`
     - attribute
     - real
     - photospheric luminosity [:math:`\erg\,\second^{-1}`]
   * - :code:`r`
     - :math:`r`
     - dataset
     - real (:code:`n`)
     - radial coordinate [:math:`\cm`]
   * - :code:`w`
     - :math:`\frac{M_{r}}{\Mstar-M_{r}}`
     - dataset
     - real (:code:`n`)
     - transformed interior mass
   * - :code:`L_r`
     - :math:`L_{r}`
     - dataset
     - real (:code:`n`)
     - interior luminosity [:math:`\erg\,\second^{-1}`]
   * - :code:`p`
     - :math:`P`
     - dataset
     - real (:code:`n`)
     - total pressure [:math:`\barye`]
   * - :code:`rho`
     - :math:`\rho`
     - dataset
     - real (:code:`n`)
     - density [:math:`\gram\,\cm^{-3}`]
   * - :code:`T`
     - :math:`T`
     - dataset
     - real (:code:`n`)
     - temperature [:math:`\kelvin`]
   * - :code:`N2`
     - :math:`N^{2}`
     - dataset
     - real (:code:`n`)
     - Brunt-Väisälä frequency squared [:math:`\second^{-2}`]
   * - :code:`Gamma_1`
     - :math:`\Gamma_{1}`
     - dataset
     - real (:code:`n`)
     - adiabatic exponent
   * - :code:`nabla_ad`
     - :math:`\nabla_{\rm ad}`
     - dataset
     - real (:code:`n`)
     - adiabatic temperature gradient
   * - :code:`delta`
     - :math:`\upsT`
     - dataset
     - real (:code:`n`)
     - thermodynamic coefficient
   * - :code:`nabla`
     - :math:`\nabla`
     - dataset
     - real (:code:`n`)
     - temperature gradient
   * - :code:`kappa`
     - :math:`\kappa`
     - dataset
     - real (:code:`n`)
     - opacity [:math:`\cm^{2}\,\gram^{-1}`]
   * - :code:`kappa_T`
     - :math:`\kapT`
     - dataset
     - real (:code:`n`)
     - opacity partial
   * - :code:`kappa_rho`
     - :math:`\kaprho`
     - dataset
     - real (:code:`n`)
     - opacity partial
   * - :code:`epsilon`
     - :math:`\epsilon`
     - dataset
     - real (:code:`n`)
     - total energy generation rate [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`epsilon_T`
     - :math:`\epsnuc\,\epsnucT`
     - dataset
     - real (:code:`n`)
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`epsilon_rho`
     - :math:`\epsnuc\,\epsnucrho`
     - dataset
     - real (:code:`n`)
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`Omega_rot`
     - :math:`\Orot`
     - dataset
     - real (:code:`n`)
     - rotation angular frequency [:math:`\radian\,\second^{-1}`]
