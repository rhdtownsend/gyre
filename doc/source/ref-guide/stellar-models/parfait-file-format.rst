.. _parfait-file-format:

PARFAIT File Format
===================

Files in PARFAIT format store HDF5 data describing a parfait stellar
model. A parfait model represents the structure of a star via a set of
concentric, spherical, uniform-density shells; see
:ads_citet:`townsend:2023b` for further details.

The PARFAIT format adheres to the following conventions:

* All data objects are attached to the root HDF5 group (`/`)
* Real values are written with type `H5T_IEEE_F64LE` when GYRE is
  compiled in double precision (the default), and type
  `H5T_IEEE_F32LE` otherwise
* Integer values are written with type `H5T_STD_I32LE`

There are a number of versions of the PARFAIT format, distinguished by
the :code:`version` attribute in the root HDF5 group:

.. toctree::

   parfait-file-format-v1.00
