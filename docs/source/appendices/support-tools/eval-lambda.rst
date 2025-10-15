.. _support-tools-eval-lambda:

eval_lambda
===========

.. program:: eval_lambda

Synopsis
--------

.. code-block:: text

   eval_lambda OUTPUT_FILE [options]

Description
-----------

The :program:`eval_lambda` tool evaluates the eigenvalue
:math:`\lambda` appearing in Laplace's tidal equations (see the
:ref:`osc-rot` section), over a grid :math:`\{q\}` of spin parameter
values. These data are written to an output file in GYRE's
:ref:`hdf-format`, with the following data items in the root group:

.. list-table::
   :widths: 10 10 10 15 55
   :header-rows: 1

   * - Item
     - Symbol
     - Object type
     - Data type
     - Definition
   * - :code:`l`
     - :math:`\ell`
     - attribute
     - integer
     - harmonic degree
   * - :code:`m`
     - :math:`m`
     - attribute
     - integer
     - azimuthal order
   * - :code:`k`
     - :math:`k`
     - attribute
     - integer
     - meridional order (see :ads_citealp:`townsend:2003a`)
   * - :code:`rossby`
     - ---
     - attribute
     - logical
     - Rossby-mode branch flag
   * - :code:`q`
     - :math:`q`
     - dataset
     - real (:math:`n`)
     - spin parameter
   * - :code:`lambda`
     - :math:`\lambda`
     - dataset
     - real (:math:`n`)
     - eigenvalue of Laplace's tidal equation


Options
-------

.. option:: -h, --help

   Print a summary of options.

.. option:: -l, --l=L

   Harmonic degree :math:`\ell`.

.. option:: -m, --m=M

   Azimuthal order :math:`m`.

.. option:: --q-min=MIN

   Minimum :math:`q` in grid.

.. option:: --q-max=MAX

   Maximum :math:`q` in grid.

.. option:: --n=N

   Number of points in grid.

.. option:: --log

   Use logarithmic grid spacing.

.. option:: --rossby

   Consider Rossby modes.
