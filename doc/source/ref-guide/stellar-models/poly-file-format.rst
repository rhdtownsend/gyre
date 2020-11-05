.. _poly-file-format:

POLY File Format
================

Files in POLY format store HDF5 data describing a :ref:`composite polytrope
model <comp-ptrope>`. This format adheres to the following conventions:

* All data objects are attached to the root HDF5 group (`/`)
* Real values are written with type `H5T_IEEE_F64LE` when GYRE is
  compiled in double precision (the default), and type
  `H5T_IEEE_F32LE` otherwise
* Integer values are written with type `H5T_STD_I32LE`

Data items in the root HDF5 group are as follows:

.. list-table::
   :widths: 8 8 8 8 68
   :header-rows: 1

   * - Data Item
     - Variable
     - Objtype
     - Datatype
     - Definition
   * - :code:`n`
     - :math:`n`
     - attribute
     - integer
     - number of grid points
   * - :code:`n_r`
     - :math:`N`
     - attribute
     - integer
     - number of regions
   * - :code:`n_poly`
     - :math:`n_{i}`
     - attribute
     - real (:code:`n_r`)
     - poltropic indices of regions
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
   * - :code:`Gamma_1`
     - :math:`\Gamma_{1}`
     - attribute
     - real
     - first adiabatic exponent
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
     - :math:`\sderiv{\theta}{z}`
     - dataset
     - real (:code:`n`)
     - Derivative of Lane-Emden variable
