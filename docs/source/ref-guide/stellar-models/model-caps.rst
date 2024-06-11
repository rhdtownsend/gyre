.. _model-caps:

Model Capabilities
==================

Which data items are included in a given stellar model dictates what
sorts of calculation can be performed on that model by the
frontends. To this end, the different model types and file formats can
be classified according to their `capabilities` (labeled using a
single letter):

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
     - no
     - yes
     - no
   * - EVOL
     - B3
     - yes
     - yes
     - no
   * - EVOL
     - FAMDL
     - no
     - yes
     - no
   * - EVOL
     - FGONG
     - no
     - yes
     - no
   * - EVOL
     - GSM
     - yes
     - yes
     - yes
   * - EVOL
     - LOSC
     - no
     - yes
     - no
   * - EVOL
     - MESA
     - yes
     - yes
     - yes
   * - EVOL
     - OSC
     - yes
     - yes
     - yes
   * - EVOL
     - WDEC
     - no
     - yes
     - no
   * - POLY
     - POLY
     - no
     - no
     - no
   * - ANAPOLY
     - ---
     - no
     - no
     - no
   * - HOM
     - ---
     - no
     - no
     - no
