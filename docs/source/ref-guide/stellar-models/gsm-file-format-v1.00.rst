Version 1.00
------------

Data items in the root HDF5 group of version-1.00 GSM-format files are as follows:

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
   * - :code:`version`
     - ---
     - attribute
     - integer
     - 100
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
   * - :code:`M_r`
     - :math:`M_r`
     - dataset
     - real (:code:`n`)
     - interior mass [:math:`\gram`]
   * - :code:`L_r`
     - :math:`L_{r}`
     - dataset
     - real (:code:`n`)
     - interior luminosity [:math:`\erg\,\second^{-1}`]
   * - :code:`P`
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
   * - :code:`kap`
     - :math:`\kappa`
     - dataset
     - real (:code:`n`)
     - opacity [:math:`\cm^{2}\,\gram^{-1}`]
   * - :code:`kap_T`
     - :math:`\kappa\,\kapT`
     - dataset
     - real (:code:`n`)
     - opacity partial [:math:`\cm^{2}\,\gram^{-1}`]
   * - :code:`kap_rho`
     - :math:`\kappa\,\kaprho`
     - dataset
     - real (:code:`n`)
     - opacity partial [:math:`\cm^{2}\,\gram^{-1}`]
   * - :code:`eps`
     - :math:`\epsilon`
     - dataset
     - real (:code:`n`)
     - total energy generation rate [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`eps_T`
     - :math:`\epsnuc\,\epsnucT`
     - dataset
     - real (:code:`n`)
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`eps_rho`
     - :math:`\epsnuc\,\epsnucrho`
     - dataset
     - real (:code:`n`)
     - nuclear energy generation rate partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - :code:`Omega_rot`
     - :math:`\Orot`
     - dataset
     - real (:code:`n`)
     - rotation angular frequency [:math:`\radian\,\second^{-1}`]
