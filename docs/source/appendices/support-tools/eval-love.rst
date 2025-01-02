.. _support-tools-eval-love:

eval_love
=========

.. program:: eval_love

The :program:`eval_love` tool evaluates the hydrostatic Love number
:math:`k_{\ell}` of an input stellar model (see
:ads_citealp:`ogilvie:2014` for the adopted definition of
:math:`k_{\ell`).

Syntax
------

:program:`eval_love` :option:`filename` :option:`model_type` :option:`file_format` :option:`l`

Parameters
----------

.. option:: filename

   Name of input stellar model file

.. option:: model_type

   Type of input stellar model (see the :nml_n:`model_type` parameter in the :ref:`model-params` section)

.. option:: file_format

   Format of input stellar model file (see the :nml_n:`file_format` parameter in the :ref:`model-params` section)

.. option:: l

   Harmonic degree :math:`\ell`

Output
------

The hydrostatic Love number is printed to standard output.

