.. _mesa-model-files:

MESA Model Files
================

MESA-format files store data describing a MESA stellar model in an
ASCII text file (note that MESA itself refers to these files as
'GYRE-format' files. To create one of these files, set the
:nml_n:`pulse_data_format` parameter of the :nml_g:`controls` namelist
group to the value :nml_v:`GYRE`).

There are a number of versions of the MESA format, distinguished by
the initial header line:

.. toctree::

   mesa-model-files-v0.01
   mesa-model-files-v0.19
   mesa-model-files-v1.00
   mesa-model-files-v1.01
