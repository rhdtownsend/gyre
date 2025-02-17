Version 0.01
------------

The first line of version-0.01 MESA-format files is a header with the following columns:

.. list-table::
   :widths: 10 10 10 70
   :header-rows: 1

   * - Column
     - Symbol
     - Datatype
     - Definition
   * - 1
     - :math:`N`
     - integer
     - number of grid points
   * - 2
     - :math:`\Mstar`
     - real
     - stellar mass [:math:`\gram`]
   * - 3
     - :math:`\Rstar`
     - real
     - photospheric radius [:math:`\cm`]
   * - 4
     - :math:`\Lstar`
     - real
     - photospheric luminosity [:math:`\erg\,\second^{-1}`]

The subsequent :math:`N` lines contain the model data, one line per
grid point extending from the center to the surface, with the
following columns:

.. list-table::
   :widths: 10 10 10 70
   :header-rows: 1

   * - Column
     - Symbol
     - Datatype
     - Definition
   * - 1
     - :math:`k`
     - integer
     - grid point index (:math:`k=1,\ldots,N`)
   * - 2
     - :math:`r`
     - real
     - radial coordinate [:math:`\cm`]
   * - 3
     - :math:`\frac{M_{r}}{\Mstar-M_{r}}`
     - real
     - transformed interior mass
   * - 4
     - :math:`L_{r}`
     - real
     - interior luminosity [:math:`\erg\,\second^{-1}`]
   * - 5
     - :math:`P`
     - real
     - total pressure [:math:`\barye`]
   * - 6
     - :math:`T`
     - real
     - temperature [:math:`\kelvin`]
   * - 7
     - :math:`\rho`
     - real
     - density [:math:`\gram\,\cm^{-3}`]
   * - 8
     - :math:`\nabla`
     - real
     - temperature gradient
   * - 9
     - :math:`N^{2}`
     - real
     - Brunt-Väisälä frequency squared [:math:`\second^{-2}`]
   * - 10
     - :math:`\cV`
     - real
     - specific heat at constant volume [:math:`\erg\,\gram^{-1}\,\kelvin^{-1}`)
   * - 11
     - :math:`\cP`
     - real
     - specific heat at constant pressure [:math:`\erg\,\gram^{-1}\,\kelvin^{-1}`)
   * - 12
     - :math:`\chi_{T}`
     - real
     - equation-of-state partial :math:`(\spderiv{\ln P}{\ln T})_{\rho}`
   * - 13
     - :math:`\chi_{\rho}`
     - real
     - equation-of-state partial :math:`(\spderiv{\ln P}{\ln \rho})_{T}`
   * - 14
     - :math:`\kappa`
     - real
     - opacity [:math:`\cm^{2}\,\gram^{-1}`]
   * - 15
     - :math:`\kapT`
     - real
     - opacity partial
   * - 16
     - :math:`\kaprho`
     - real
     - opacity partial
   * - 17
     - :math:`\epsilon`
     - real
     - total energy generation rate [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - 18
     - :math:`\epsnuc\,\epsnucT`
     - real
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - 19
     - :math:`\epsnuc\,\epsnucrho`
     - real
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
