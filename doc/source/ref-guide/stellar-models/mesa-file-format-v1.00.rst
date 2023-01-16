Version 1.00
------------

The first line of version-1.00 MESA-format files is a header with the following columns:

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
     - :math:`M`
     - real
     - stellar mass [:math:`\gram`]
   * - 3
     - :math:`R`
     - real
     - photospheric radius [:math:`\cm`]
   * - 4
     - :math:`L`
     - real
     - photospheric luminosity [:math:`\erg\,\second^{-1}`]
   * - 5
     - 100
     - integer
     - version number

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
     - :math:`M_{r}`
     - real
     - interior mass [:math:`\gram`]
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
     - dimensionless temperature gradient
   * - 9
     - :math:`N^{2}`
     - real
     - Brunt-Väisälä frequency squared [:math:`\second^{-2}`]
   * - 10
     - :math:`\Gamma_{1}`
     - real
     - adiabatic exponent
   * - 11
     - :math:`\nabla_{\rm ad}`
     - real
     - adiabatic temperature gradient
   * - 12
     - :math:`\upsT`
     - real
     - thermodynamic coefficient
   * - 13
     - :math:`\kappa`
     - real
     - opacity [:math:`\cm^{2}\,\gram^{-1}`]
   * - 14
     - :math:`\kappa\,\kapT`
     - real
     - opacity partial [:math:`\cm^{2}\,\gram^{-1}`]
   * - 15
     - :math:`\kappa\,\kaprho`
     - real
     - opacity partial [:math:`\cm^{2}\,\gram^{-1}`]
   * - 16
     - :math:`\epsilon`
     - real
     - total energy generation rate [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - 17
     - :math:`\epsnuc\,\epsnucT`
     - real
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - 18
     - :math:`\epsnuc\,\epsnucrho`
     - real
     - nuclear energy generation partial [:math:`\erg\,\second^{-1}\,\gram^{-1}`]
   * - 19
     - :math:`\Omega_{\rm rot}`
     - real
     - rotation angular velocity [:math:`\radian\,\second^{-1}`]
