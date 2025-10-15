.. _support-tools-poly-to-fgong:

poly_to_fgong
=============

.. program:: poly_to_fgong

Synopsis
--------

.. code-block:: text

   poly_to_fgong POLY_FILE FGONG_FILE [options]

Description
-----------

The :program:`poly_to_fgong` tool converts :ref:`polytropic stellar
models <poly-models>` in the :ref:`POLY <poly-file-format>` format, to
:ref:`evolutionary models <evol-models>` in the FGONG format. In this
conversion, a stellar mass :math:`M = 1\,\Msun` and radius :math:`R =
1\,\Rsun` are presumed.

Options
-------

.. option:: -h, --help

   Print a summary of options.

.. option:: --drop_outer

   Drop the outer-most point (which is singular when the surface
   pressure/density vanish).
