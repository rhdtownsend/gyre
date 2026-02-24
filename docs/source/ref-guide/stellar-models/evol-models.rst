.. _evol-models:

Evolutionary Models
===================

.. nml:group:: model
   :no-target:

Setting the :nml:option:`model_type` option of the :nml:group:`model`
namelist group to :nml:value:`'EVOL'` tells the frontend to read the
equilibrium stellar model from a file created by a stellar evolution
code (e.g., `MESA <mesa_>`__).

Supported Formats
-----------------

The format of the model file is specified by the
:nml:option:`file_format` option of the :nml:group:`model` namelist
group (see the :ref:`model-group` section). Possible choices are
summarized in the table below.

.. list-table::
   :widths: 20 80
   :header-rows: 1
   :align: center

   * - :nml:option:`file_format`
     - Description
   * - :nml:value:`'AMDL'`
     - Binary file describing an evolutionary model in AMDL format,
       as reverse engineered from the ADIPLS stellar oscillation code
       :ads_citep:`christensen-dalsgaard:2008`
   * - :nml:value:`'B3'`
     - HDF5 file describing an evolutionary model in B3 format. This
       format is for testing purposes only, and will eventually be
       superseded and/or removed
   * - :nml:value:`'FAMDL'`
     - Text file describing an evolutionary model in FAMDL format, as
       specified in the :download:`CoRoT/ESTA File Formats
       <corot-esta-file-formats.pdf>` document
   * - :nml:value:`'FGONG'`
     - Text file describing an evolutionary model in FGONG format, as
       specified in the updated :download:`FGONG Format
       <fgong-file-format.pdf>` document
   * - :nml:value:`'GSM'`
     - HDF5 file describing an evolutionary model in GYRE Stellar
       Model (GSM) format, as specified in the :ref:`gsm-file-format`
       section
   * - :nml:value:`'MESA'`
     - Text file describing an evolutionary model in MESA format, as
       specified in the :ref:`mesa-file-format` section
   * - :nml:value:`'LOSC'`
     - Text file describing an evolutionary model in the revised LOSC
       format
   * - :nml:value:`'OSC'`
     - Text file describing an evolutionary model in OSC format, as
       specified in the :download:`CoRoT/ESTA File Formats
       <corot-esta-file-formats.pdf>` document)
   * - :nml:value:`'WDEC'`
     - Text file describing an evolutionary model in WDEC format, as
       specified in :ads_citet:`bischoff-kim:2018`

Interpolation
-------------

Cubic spline interpolation is used to evaluate data between model grid
points. The :nml:option:`deriv_type` option controls how the spline
derivatives are set up.

.. _evol-models-double:

Double Points
-------------

If a model contains a pair of adjacent points with the same radial
coordinate :math:`r`, this pair is treated as a double point
representing a discontinuity in the density and some other
thermodynamic quantities (but not the pressure or temperature). GYRE
does not attempt to interpolate across double points, but instead
handles them properly when solving equations through the use of
:ref:`internal boundary conditions <osc-bound-conds>`.
