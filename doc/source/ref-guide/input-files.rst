.. _namelist-input-files:

********************
Namelist Input Files
********************

GYRE executables read parameters from input files that defines a
number of Fortran namelist groups. The table below lists the groups,
and for each exeutable indicates how many of these groups can appear
in a valid input file.

.. toctree::
   :maxdepth: 2
   :hidden:

   input-files/constants
   input-files/grid-params
   input-files/model-params
   input-files/mode-params
   input-files/num-params
   input-files/orbit-params
   input-files/osc-params
   input-files/output-params
   input-files/rot-params
   input-files/scan-params
   input-files/tide-params

.. list-table::
   :header-rows: 1
   :widths: 30 30 20 20

   * - Description
     - Namelist group name
     - For gyre
     - For gyre_tides
   * - :ref:`constants`
     - :nml_g:`constants`
     - 1
     - 1
   * - :ref:`grid-params`
     - :nml_g:`grid`
     - :math:`\geq 1`\ [#last]_
     - :math:`\geq 1`\ [#last]_
   * - :ref:`mode-params`
     - :nml_g:`mode`
     - :math:`\geq 1`
     - ---
   * - :ref:`model-params`
     - :nml_g:`model`
     - 1
     - 1
   * - :ref:`num-params`
     - :nml_g:`num`
     - :math:`\geq 1`\ [#last]_
     - :math:`\geq 1`\ [#last]_
   * - :ref:`orbit-params`
     - :nml_g:`orbit`
     - ---
     - :math:`\geq 1`\ [#last]_
   * - :ref:`osc-params`
     - :nml_g:`osc`
     - :math:`\geq 1`\ [#last]_
     - :math:`\geq 1`\ [#last]_
   * - :ref:`output-params`
     - :nml_g:`ad_output`
     - 1
     - ---
   * - 
     - :nml_g:`nad_output`
     - 1
     - ---
   * - 
     - :nml_g:`tides_output`
     - ---
     - 1
   * - :ref:`rot-params`
     - :nml_g:`rot`
     - :math:`\geq 1`\ [#last]_
     - :math:`\geq 1`\ [#last]_
   * - :ref:`scan-params`
     - :nml_g:`scan`
     - :math:`\geq 1`
     - ---
   * - :ref:`tide-params`
     - :nml_g:`tide`
     - ---
     - :math:`\geq 1`
       
.. rubric:: Footnotes

.. [#last] While the input file can contain one or more of the indicated namelist group, only the last (tag-matching) one is used.
