.. _mesa-file-format:

MESA File Format
================

Files in MESA format store ASCII text data describing a MESA stellar
model. Note that MESA itself refers to these files as 'GYRE-format'
files; to create one of these files in MESA, set the
:nml_n:`pulse_data_format` parameter of the :nml_g:`controls` namelist
group to the value :nml_v:`GYRE`.

There are a number of versions of the MESA format, distinguished by
the initial header line:

.. toctree::

   mesa-file-format-v0.01
   mesa-file-format-v0.19
   mesa-file-format-v1.00
   mesa-file-format-v1.01
