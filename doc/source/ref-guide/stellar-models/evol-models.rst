.. _evol-models:

Evolutionary Models
===================

Supported Formats
-----------------

Evolutionary models are read from a file created by a separate stellar
evolution code. The format of this file is specified by the
:nml_n:`file_format` parameter of the :nml_g:`model` namelist group
(see the :ref:`model-params` section). Possible choices are summarized
in the table below.

.. list-table::
   :widths: 20 80
   :header-rows: 1
   :align: center

   * - :nml_n:`file_format`
     - Description
   * - :nml_v:`'AMDL'`
     - Binary file describing an evolutionary model in AMDL format,
       as reverse engineered from the ADIPLS stellar oscillation code
       :ads_citep:`christensen-dalsgaard:2008`
   * - :nml_v:`'B3'`
     - HDF5 file describing an evolutionary model in B3 format. This
       format is for testing purposes only, and will eventually be
       superseded and/or removed
   * - :nml_v:`'FAMDL'`
     - Text file describing an evolutionary model in FAMDL format, as
       specified in the :download:`CoRoT/ESTA File Formats
       <corot-esta-file-formats.pdf>` document
   * - :nml_v:`'FGONG'`
     - Text file describing an evolutionary model in FGONG format, as
       specified in the updated :download:`FGONG Format
       <fgong-file-format.pdf>` document
   * - :nml_v:`'GSM'`
     - HDF5 file describing an evolutionary model in GYRE Stellar
       Model (GSM) format, as specified in the :ref:`gsm-file-format`
       section
   * - :nml_v:`'MESA'`
     - Text file describing an evolutionary model in MESA format, as
       specified in the :ref:`mesa-file-format` section
   * - :nml_v:`'LOSC'`
     - Text file describing an evolutionary model in the revised LOSC
       format
   * - :nml_v:`'OSC'`
     - Text file describing an evolutionary model in OSC format, as
       specified in the :download:`CoRoT/ESTA File Formats
       <corot-esta-file-formats.pdf>` document)
   * - :nml_v:`'WDEC'`
     - Text file describing an evolutionary model in WDEC format, as
       specified in :ads_citet:`bischoff-kim:2018`

Interpolation
-------------
  
Cubic spline interpolation is used to evaluate data between model grid
points. The :nml_n:`deriv_type` parameter in the :nml_g:`model`
namelist group controls how the spline derivatives are set up.

.. _evol-models-double:

Double Points
-------------

If a model contains a pair of adjacent points with the same radial
coordinate :math:`r`, this pair is treated as a double point
representing a discontinuity in the density and some other
thermodynamic quantities (but not the pressure or temperature). GYRE
does not attempt to interpolate across double points, but does handle
them properly when solving the oscillation equations through the use
of :ref:`jump conditions <dimless-form-jump>`.
