.. _evol-models:

Evolutionary Models
===================

Evolutionary models are read from a file created by a separate stellar
evolution code. The format of this file is specified by the
:nml_n:`file_format` parameter of the :nml_g:`model` namelist group
(see the :ref:`model-params` section). Possible choices are:

:nml_v:`AMDL`
  Binary file describing an evolutionary model in ASTEC format (as
  reverse engineered from the ADIPLS stellar oscillation code; see
  :ads_citealp:`christensen-dalsgaard:2008`).

:nml_v:`FAMDL`
  Text file describing an evolutionary model in FAMDL format (as
  specified in the :download:`CoRoT/ESTA File Formats
  <corot-esta-file-formats.pdf>` document).

:nml_v:`FGONG`
  Text file describing an evolutionary model in FGONG format (as
  specified in the updated :download:`FGONG Format
  <fgong-file-format.pdf>` document).

:nml_v:`GSM`
  HDF5 file describing an evolutionary model in GYRE
  Stellar Model (GSM) format (as specified in the :ref:`gsm-file-format` section).

:nml_v:`MESA`
  Text file describing an evolutionary model in MESA format (as
  specified in the :ref:`mesa-file-format` section).
  
:nml_v:`B3`
  HDF5 file describing an evolutionary model in B3 format. This format
  is for testing purposes only, and will eventually be superseded and/or
  removed.

:nml_v:`LOSC`
  Text file describing an evolutionary model in the revised LOSC
  format.

:nml_v:`OSC`
  Text file describing an evolutionary model in OSC format (as
  specified in the :download:`CoRoT/ESTA File Formats
  <corot-esta-file-formats.pdf>` document).

:nml_v:`WDEC`
  Text file describing an evolutionary model in WDEC format (see
  :ads_citealp:`bischoff-kim:2018`)
  
For all of these model formats, cubic spline interpolation is used to
evaluate data between model grid points. The :nml_n:`deriv_type`
parameter in the :nml_g:`model` namelist group controls how the spline
derivatives are set up.
