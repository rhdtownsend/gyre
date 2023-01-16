.. _gsm-file-format:

GSM File Format
===============

Files in GSM (GYRE Stellar Model) format store HDF5 data describing a
stellar model. This format is intended as a portable,
storage-efficient alternative to the :ref:`mesa-file-format`. To
create one of these files in MESA (from revision 21.12.1 onward), set
the :nml_n:`pulse_data_format` parameter of the :nml_g:`controls`
namelist group to the value :nml_v:`'GSM'`.

The GSM format adheres to the following conventions:

* All data objects are attached to the root HDF5 group (`/`)
* Real values are written with type `H5T_IEEE_F64LE` when GYRE is
  compiled in double precision (the default), and type
  `H5T_IEEE_F32LE` otherwise
* Integer values are written with type `H5T_STD_I32LE`

There are a number of versions of the GSM format, distinguished by the
:code:`version` attribute in the root HDF5 group:

.. toctree::

   gsm-file-format-v0.00
   gsm-file-format-v1.00
   gsm-file-format-v1.10
