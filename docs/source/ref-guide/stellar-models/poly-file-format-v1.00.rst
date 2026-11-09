Version 1.00
------------

The data objects appearing in the root HDF5 group of version-1.00
POLY-format files are as follows:

.. list-table::
   :widths: 10 10 15 20 45
   :header-rows: 1

   * - Object name
     - Variable
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
   * - :code:`n_r`
     - :math:`\nreg`
     - attribute
     - integer
     - number of regions
   * - :code:`n_poly`
     - :math:`n_{i}`
     - attribute
     - real (:code:`n_r`)
     - polytropic indices of regions
   * - :code:`z_b`
     - :math:`z_{i-1/2}`
     - attribute
     - real (:code:`n_r-1`)
     - radial coordinates of region boundaries
   * - :code:`Delta_b`
     - :math:`\Delta_{i-1/2}`
     - attribute
     - real (:code:`n_r-1`)
     - log density jump of region boundaries
   * - :code:`z`
     - :math:`z`
     - dataset
     - real (:code:`n`)
     - polytropic radial coordinate
   * - :code:`theta`
     - :math:`\theta`
     - dataset
     - real (:code:`n`)
     - Lane-Emden variable
   * - :code:`dtheta`
     - :math:`\theta'`
     - dataset
     - real (:code:`n`)
     - Derivative of Lane-Emden variable
   * - :code:`Gamma_1`
     - :math:`\Gamma_{1}`
     - dataset
     - real (:code:`n`)
     - first adiabatic exponent
   * - :code:`Gamma_A`
     - :math:`\Gamma_{A}`
     - dataset
     - real (:code:`n`)
     - first adiabatic exponent used when evaluating :math:`\As` structure coefficient
       
