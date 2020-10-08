Version 0.01
------------

The first line of version-0.01 MESA-format files is a header with the following columns:

.. list-table::
   :widths: 15 15 15 55
   :header-rows: 1

   * - Column
     - Variable
     - Datatype
     - Definition
   * - 1
     - :math:`n`
     - integer
     - number of grid points
   * - 2
     - :math:`M`
     - real
     - stellar mass (:math:`\gram`)
   * - 3
     - :math:`R`
     - real
     - photospheric radius (:math:`\cm`)
   * - 4
     - :math:`L`
     - real
     - photospheric luminosity (:math:`\erg\,\second^{-1}`)

The subsequent :math:`n` lines contain the model data, one line per
grid point extending from the center to the surface, with the
following columns:

.. list-table::
   :widths: 10 10 10 70
   :header-rows: 1

   * - Column
     - Variable
     - Datatype
     - Definition
   * - 1
     - :math:`k`
     - integer
     - grid point index (:math:`k=1,\ldots,n`)
   * - 2
     - :math:`r`
     - real
     - radial coordinate (:math:`\cm`)
   * - 3
     - :math:`w`
     - real
     - transformed mass coordinate :math:`M_{r}/(M-M_{r})`
   * - 4
     - :math:`L_{r}`
     - real
     - luminosity (:math:`\erg\,\second^{-1}`)
   * - 5
     - :math:`P`
     - real
     - total pressure (:math:`\dyne\,\cm^{-2}`)
   * - 6
     - :math:`T`
     - real
     - temperature (:math:`\kelvin`)
   * - 7
     - :math:`\rho`
     - real
     - density (:math:`\gram\,\cm^{-3}`)
   * - 8
     - :math:`\nabla`
     - real
     - dimensionless temperature gradient :math:`\sderiv{\ln T}{\ln P}`
   * - 9
     - :math:`N^{2}`
     - real
     - Brunt-Väisälä frequency squared (:math:`\second^{-2}`)
   * - 10
     - :math:`c_{v}`
     - real
     - specific heat at constant volume (:math:`\erg\,\gram^{-1}\,\kelvin^{-1}`)
   * - 11
     - :math:`c_{p}`
     - real
     - specific heat at constant pressure (:math:`\erg\,\gram^{-1}\,\kelvin^{-1}`)
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
     - opacity (:math:`\cm^{2}\,\gram^{-1}`)
   * - 15
     - :math:`\kapT`
     - real
     - opacity partial :math:`(\spderiv{\ln \kappa}{\ln T})_{\rho}`
   * - 16
     - :math:`\kaprho`
     - real
     - opacity partial :math:`(\spderiv{\ln \kappa}{\ln \rho})_{T}`
   * - 17
     - :math:`\epsilon`
     - real
     - total energy generation rate (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - 18
     - :math:`\epsnuc\,\epsT`
     - real
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln T})_{\rho}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - 19
     - :math:`\epsnuc\,\epsrho`
     - real
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln \rho})_{T}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
