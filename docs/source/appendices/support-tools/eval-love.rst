.. _support-tools-eval-love:

eval_love
=========

.. program:: eval_love

Synopsis
--------

.. code-block:: text

   eval_love MODEL_FILE [options]

Description
-----------

The :program:`eval_love` tool evaluates the hydrostatic Love number
:math:`k_{\ell}` of an input stellar model (see
:ads_citealp:`ogilvie:2014` for the adopted definition of
:math:`k_{\ell}`), and writes it to standard output.

Options
-------

.. option:: -h, --help

   Print a summary of options.

.. option:: --model-type=TYPE

   Type of input stellar model (see the :nml:option:`model_type <model.model_type>` option)

.. option:: --file-format=FORMAT

   Format of input stellar model file (see the :nml:option:`file_format <model.file_format>` option)

.. option:: -l, --l=L

   Harmonic degree :math:`\ell`.
