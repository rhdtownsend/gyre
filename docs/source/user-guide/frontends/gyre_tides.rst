.. _frontends-gyre_tides:

gyre_tides
==========

.. image:: gyre_tides-flow.drawio.svg
   :align: right
   :alt: Flow of execution in the gyre_tides frontend.

The :program:`gyre_tides` frontend calculates the response of
a stellar model to tidal forcing by a orbiting point-mass
companion. The general flow of execution is outlined in the chart to
the right. After reading the :ref:`namelist input file
<namelist-input-files>` and the :ref:`model <stellar-models>`,
:program:`gyre_tides` loops over :nml:group:`tide` namelist groups,
processing each in turn.

For a given group, :program:`gyre_tides` solves for the response of
the star to the superposition of partial tidal potentials
:math:`\PhiTlmk` (see the :ref:`osc-tidal` section). The response
wavefunctions and other data associated with an individual partial
potential are optionally written to a :ref:`detail file
<detail-files>`.  At the end of the run, response data from all
partial potentials (across all :nml:group:`tide` groups) are optionally
written to a :ref:`summary file <summary-files>`.

The table below lists which namelist groups, and in what number,
should appear in namelist input files for :program:`gyre_tides`.

.. list-table::
   :header-rows: 1
   :widths: 30 30 20

   * - Description
     - Namelist group name
     - Number
   * - :ref:`constants-group`
     - :nml:group:`constants`
     - 1
   * - :ref:`grid-group`
     - :nml:group:`grid`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`model-group`
     - :nml:group:`model`
     - 1
   * - :ref:`num-group`
     - :nml:group:`num`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`orbit-group`
     - :nml:group:`orbit`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`osc-group`
     - :nml:group:`osc`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`output-groups`
     - :nml:group:`tide_output`
     - 1
   * - :ref:`rot-group`
     - :nml:group:`rot`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`tide-group`
     - :nml:group:`tide`
     - :math:`\geq 1`

.. rubric:: Footnotes

.. [#last] While the input file can contain one or more of the
           indicated namelist group, only the last (:ref:`tag-matching
           <working-with-tags>`) one is used.
