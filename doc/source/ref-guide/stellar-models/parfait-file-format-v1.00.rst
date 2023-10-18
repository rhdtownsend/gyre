Version 1.00
------------

Data items in the root HDF5 group of version-1.00 PARFAIT-format files are as follows:

.. list-table::
   :widths: 5 10 10 15 60
   :header-rows: 1

   * - Item
     - Symbol
     - Object type
     - Data type
     - Definition
   * - :code:`type`
     - ---
     - attribute
     - string
     - :code:`'PARFAIT'`
   * - :code:`version`
     - ---
     - attribute
     - integer
     - 100
   * - :code:`N`
     - :math:`N`
     - attribute
     - integer
     - number of shells
   * - :code:`y_c`
     - :math:`m_{1/2}\,/\,M`
     - attribute
     - real
     - normalized central mass coordinate
   * - :code:`z_s`
     - :math:`P_{N+1/2}\,/\,(GM^{2}/R^{4})`
     - attribute
     - real
     - normalized surface pressure
   * - :code:`x`
     - :math:`r/R`
     - dataset
     - real (:code:`N+1`)
     - normalized shell interface radial coordinate
   * - :code:`d`
     - :math:`\rho\,/\,(M/R^{3})`
     - dataset
     - real (:code:`N`)
     - normalized shell density
   * - :code:`Gamma_1`
     - :math:`\Gamma_{1}`
     - dataset
     - real (:code:`N`)
     - shell adiabatic exponent
