.. _model-caps:

Model Capabilities
==================

Which data items are included in a given stellar model dictates what
sorts of GYRE calculation can be performed on that model. To this end,
the different model types and file formats can be classified according
to their `capabilities` (labeled using a single letter):

N
  The model supports non-adiabatic calculations.

D
  The model supports dimensioned results.

R
  The model supports differential rotational.

The table below summarizes the different capabilities of each
model-type and file-format combination.

.. list-table::
   :widths: 30 30 10 10 10
   :header-rows: 1
   :align: center

   * - Model Type
     - File Format
     - N
     - D
     - R
   * - EVOL
     - AMDL
     - 
     - X
     - 
   * - EVOL
     - B3
     - X
     - X
     -
   * - EVOL
     - FAMDL
     -
     - X
     - 
   * - EVOL
     - FGONG
     -
     - X
     -
   * - EVOL
     - GSM
     - X
     - X
     - X
   * - EVOL
     - MESA
     - X
     - X
     - X
   * - EVOL
     - LOSC
     -
     - X
     - 
   * - EVOL
     - OSC
     - X
     - X
     - X
   * - EVOL
     - WDEC
     -
     - X
     - 
   * - POLY
     - POLY
     -
     -
     -
   * - HOM
     - ---
     -
     -
     -
