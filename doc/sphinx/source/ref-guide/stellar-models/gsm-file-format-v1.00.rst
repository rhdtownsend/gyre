Version 1.00
------------

Data items in the root HDF5 group of version-1.00 GSM-format files are as follows:

.. list-table::
   :widths: 8 8 8 66
   :header-rows: 1

   * - Data Item
     - Variable
     - Datatype
     - Definition
   * - :code:`n`
     - :math:`n`
     - integer scalar
     - number of grid points
   * - :code:`version`
     - :math:`\text{version} \times 100`
     - integer scalar
     - 100
   * - :code:`M_star`
     - :math:`M`
     - real scalar
     - stellar mass (:math:`\gram`)
   * - :code:`R_star`
     - :math:`R`
     - real scalar
     - photospheric radius (:math:`\cm`)
   * - :code:`L_star`
     - :math:`L`
     - real scalar
     - photospheric luminosity (:math:`\erg\,\second^{-1}`)
   * - :code:`r`
     - :math:`r`
     - real array
     - radial coordinate (:math:`\cm`)
   * - :code:`M_r`
     - :math:`M_r`
     - real array
     - interior mass (:math:`\gram`)
   * - :code:`L_r`
     - :math:`L_{r}`
     - real array
     - interior luminosity (:math:`\erg\,\second^{-1}`)
   * - :code:`P`
     - :math:`P`
     - real array
     - total pressure (:math:`\barye`)
   * - :code:`rho`
     - :math:`\rho`
     - real array
     - density (:math:`\gram\,\cm^{-3}`)
   * - :code:`T`
     - :math:`T`
     - real array
     - temperature (:math:`\kelvin`)
   * - :code:`N2`
     - :math:`N^{2}`
     - real array
     - Brunt-Väisälä frequency squared (:math:`\second^{-2}`)
   * - :code:`Gamma_1`
     - :math:`\Gamma_{1}`
     - real array
     - first adiabatic exponent :math:`(\spderiv{\ln P}{\ln \rho})_{\rm ad}`
   * - :code:`nabla_ad`
     - :math:`\nabla_{\rm ad}`
     - real array
     - adiabatic temperature gradient :math:`(\spderiv{\ln T}{\ln P})_{\rm ad}`
   * - :code:`delta`
     - :math:`\delta`
     - real array
     - dimensionless thermal expansion coefficient :math:`-(\spderiv{\ln \rho}{\ln T})_{P}`
   * - :code:`nabla`
     - :math:`\nabla`
     - real array
     - dimensionless temperature gradient :math:`\sderiv{\ln T}{\ln P}`
   * - :code:`kap`
     - :math:`\kappa`
     - real array
     - opacity (:math:`\cm^{2}\,\gram^{-1}`)
   * - :code:`kap_T`
     - :math:`\kappa\,\kapT`
     - real array
     - opacity partial :math:`\kappa (\spderiv{\ln \kappa}{\ln T})_{\rho}` (:math:`\cm^{2}\,\gram^{-1}`)
   * - :code:`kap_rho`
     - :math:`\kappa\,\kaprho`
     - real array
     - opacity partial :math:`\kappa (\spderiv{\ln \kappa}{\ln \rho})_{T}` (:math:`\cm^{2}\,\gram^{-1}`)
   * - :code:`eps`
     - :math:`\epsilon`
     - real array
     - total energy generation rate (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - :code:`eps_T`
     - :math:`\epsnuc\,\epsT`
     - real array
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln T})_{\rho}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - :code:`eps_rho`
     - :math:`\epsnuc\,\epsrho`
     - real array
     - nuclear energy generation rate partial :math:`\epsnuc (\spderiv{\ln \epsnuc}{\ln \rho})_{T}` (:math:`\erg\,s^{-1}\,\gram^{-1}`)
   * - :code:`Omega_rot`
     - :math:`\Omega_{\rm rot}`
     - real array
     - rotation angular velocity (:math:`\radian\,\second^{-1}`)
