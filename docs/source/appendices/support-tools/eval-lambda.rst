.. _support-tools-eval-lambda:

eval_lambda
===========

.. program:: eval_lambda

The :program:`eval_lambda` tool evaluates the eigenvalue
:math:`\lambda` appearing in Laplace's tidal equations (see the
:ref:`osc-rot` section), over a grid :math:`\{q\}` of spin parameter
values.

Syntax
------

:program:`eval_lambda` :option:`l` :option:`m` :option:`q_min` :option:`q_max` :option:`n_q` :option:`log_q` :option:`rossby` :option:`filename`

Parameters
----------

.. option:: l

   Harmonic degree :math:`\ell`

.. option:: m

   Azimuthal order :math:`m`.

.. option:: q_min

   Start value in spin parameter grid

.. option:: q_max

   End value in spin parameter grid

.. option:: n_q

   Number of points in spin parameter grid

.. option:: log_q

   Flag to use logarithmic spacing in spin parameter grid

.. option:: rossby

   Flag to consider the Rossby-mode branch of Laplace's tidal equations

.. option:: filename

   Name of output file (see below)

Output
------

The output file created by :program:`eval_lambda` is in GYRE's
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
     - real (:option:`n_q`)
     - spin parameter
   * - :code:`lambda`
     - :math:`\lambda`
     - dataset
     - real (:option:`n_q`)
     - eigenvalue of Laplace's tidal equation
