.. _support-tools-poly-to-fgong:

poly_to_fgong
=============

.. program:: poly_to_fgong

The :program:`poly_to_fgong` tool converts :ref:`polytropic stellar
models <poly-models>` in the :ref:`POLY <poly-file-format>` format, to
:ref:`evolutionary models <evol-models>` in the FGONG format. In this
conversion, a stellar mass :math:`M = 1\,\Msun` and radius :math:`R =
1\,\Rsun` are presumed.

Syntax
------

:program:`poly_to_fgong` :option:`in_filename` :option:`out_filename` :option:`drop_outer`

Parameters
----------

.. option:: in_filename

   File name of input POLY file

.. option:: out_filename

   File name of output FGONG file

.. option:: drop_outer

   Flag to drop outer-most point (which is singular when the surface pressure/density vanish)

Output
------

The output is written to the :option:`out_filename` file.

