Version 0.19
------------

The first line of version-0.19 MESA-format files is a header with the following columns:

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
   * - 5
     - 19
     - integer
     - version number :math:`\times 100`

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
     - :math:`\Gamma_{1}`
     - real
     - first adiabatic exponent :math:`(\spderiv{\ln P}{\ln \rho})_{\rm ad}`
   * - 11
     - :math:`\nabla_{\rm ad}`
     - real
     - adiabatic temperature gradient :math:`(\spderiv{\ln T}{\ln P})_{\rm ad}`
   * - 12
     - :math:`\delta`
     - real
     - dimensionless thermal expansion coefficient :math:`-(\spderiv{\ln \rho}{\ln T})_{P}`
   * - 13
     - :math:`\kappa`
     - real
     - opacity (:math:`\cm^{2}\,\gram^{-1}`)
   * - 14
     - :math:`\kapT`
     - real
     - opacity partial :math:`(\spderiv{\ln \kappa}{\ln T})_{\rho}`
   * - 15
     - :math:`\kaprho`
     - real
     - opacity partial :math:`(\spderiv{\ln \kappa}{\ln \rho})_{T}`
   * - 16
     - :math:`\epsilon`
     - real
     - total energy generation rate (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - 17
     - :math:`\epsnuc\,\epsT`
     - real
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln T})_{\rho}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - 18
     - :math:`\epsnuc\,\epsrho`
     - real
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln \rho})_{T}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - 19
     - :math:`\Omega_{\rm rot}`
     - real
     - rotation angular velocity (:math:`\radian\,\second^{-1}`)
