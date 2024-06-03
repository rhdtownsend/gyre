.. _eval-lambda:

****************************
Evaluating Tidal Eigenvalues
****************************

This appendix describes the :program:`eval_lambda` executable, which
evaluates the eigenvalue :math:`\lambda` appearing in Laplace's tidal
equations (see the :ref:`osc-rot` section). This executable is
used for the calculations presented in :ads_citet:`townsend:2020`.

Installation
============

:program:`eval_lambda` is automatically compiled when GYRE is built,
and installed in the :file:`{$GYRE_DIR}/bin` directory (see the main
:ref:`installation` chapter).

Running
=======

Unlike most other GYRE executables, the parameters for
:program:`eval_lambda` are supplied directly on the command line, with
the syntax

.. prompt:: bash

   ./eval_lambda l m q_min q_max n_q log_q rossby filename

This evaluates :math:`\lambda` for harmonic degree :math:`\ell` and
azimuthal order :math:`m` on a grid
:math:`\{q_{1},q_{2},\ldots,q_{N}\}` in the spin parameter, writing
the results to the file :file:`filename`. If the flag :code:`log_q`
has the value :code:`T` then the grid is logarithmically spaced:

.. math::

   q_{i} = 10^{(1 - w_{i}) \log q_{\rm min} + w_{i} \log q_{\rm max}},

where

.. math::

   w_{i} \equiv \frac{i-1}{N-1}.

Alternatively, if :code:`log_q` has the value :code:`F`, then the grid
is linearly spaced:

.. math::

   q_{i} = (1 - w_{i}) q_{\rm min} + w_{i} q_{\rm max}.

As a special case, when :math:`n_{q}=1`, :math:`q_{\rm min}` and
:math:`q_{\rm max}` must match, and the single :math:`q` point is set
to equal them.

If the flag :code:`rossby` has the value :code:`T`, then the Rossby
:math:`\lambda` family is evaluated; otherwise, the gravito-acoustic
family is evaluated.

The table below summarizes the mapping between the user-definable
controls appearing in the expressions above, and the corresponding
command-line parameters:

.. list-table::
   :widths: 30 30 
   :header-rows: 1

   * - Symbol
     - Parameter
   * - :math:`\ell`
     - :code:`l`
   * - :math:`m`
     - :code:`m`
   * - :math:`q_{\rm min}`
     - :code:`q_min`
   * - :math:`q_{\rm max}`
     - :code:`q_max`
   * - :math:`N`
     - :code:`n_q`
   
Interpreting Output
===================

The output file created by :program:`eval_lambda` is in GYRE's
:ref:`hdf-format`, with the following data:

:code:`l` (integer scalar)
  Harmonic degree :math:`\ell`

:code:`k` (integer scalar)
  Meridional order :math:`k` (see :ads_citealp:`townsend:2003a`)

:code:`m` (integer scalar)
  Azimuthal order :math:`m`

:code:`rossby` (logical scalar)
  Rossby family flag

:code:`q` (real array)
  Spin parameter :math:`q`

:code:`lambda` (real array)
  Eigenvalue :math:`\lambda`
