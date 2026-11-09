.. _poly-file-format:

POLY File Format
================

Files in POLY format store HDF5 data describing a :ref:`composite polytrope
model <comp-ptrope>`. This format adheres to the following conventions:

* All data objects are attached to the root HDF5 group (`/`)
* Real values are written with type `H5T_IEEE_F64LE`
* Integer values are written with type `H5T_STD_I32LE`

There are a number of versions of the GSM format, distinguished by the
:code:`version` attribute in the root group:

.. toctree::

   poly-file-format-v0.00
   poly-file-format-v1.00
