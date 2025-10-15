.. _frontends-gyre:

gyre
====

.. image:: gyre-flow.drawio.svg
   :align: right
   :alt: Flow of execution in the gyre frontend.

The :program:`gyre` frontend calculates the free-oscillation
modes of a stellar model. The general flow of execution is outlined in
the chart to the right. After reading the :ref:`namelist input file
<namelist-input-files>` and the :ref:`model <stellar-models>`,
:program:`gyre` loops over :nml_g:`mode` namelist groups,
processing each in turn.

For a given group, :program:`gyre` searches over a range of
oscillation frequencies for modes with a specific harmonic degree
:math:`\ell` and azimuthal order :math:`m`. With each mode found, the
eigenfrequency, eigenfunctions and other data are optionally written
to a :ref:`detail file <detail-files>`.  At the end of the run,
response data from all modes found (across all :nml_g:`mode` groups)
are optionally written to a :ref:`summary file <summary-files>`.

The table below lists which namelist groups, and in what number,
should appear in namelist input files for :program:`gyre`.

.. list-table::
   :header-rows: 1
   :widths: 30 30 20

   * - Description
     - Namelist group name
     - Count
   * - :ref:`constants`
     - :nml_g:`constants`
     - 1
   * - :ref:`grid-params`
     - :nml_g:`grid`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`mode-params`
     - :nml_g:`mode`
     - :math:`\geq 1`
   * - :ref:`model-params`
     - :nml_g:`model`
     - 1
   * - :ref:`num-params`
     - :nml_g:`num`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`osc-params`
     - :nml_g:`osc`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`output-params`
     - :nml_g:`ad_output`
     - 1
   * -
     - :nml_g:`nad_output`
     - 1
   * - :ref:`rot-params`
     - :nml_g:`rot`
     - :math:`\geq 1`\ [#last]_
   * - :ref:`scan-params`
     - :nml_g:`scan`
     - :math:`\geq 1`

.. rubric:: Footnotes

.. [#last] While the input file can contain one or more of the
           indicated namelist group, only the last (:ref:`tag-matching
           <working-with-tags>`) one is used.
