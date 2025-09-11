.. _support-tools-poly-to-txt:

poly_to_txt
===========

.. program:: poly_to_txt

Synopsis
--------

.. code-block:: text

   poly_to_txt POLY_FILE TXT_FILE [options]

Description
-----------

The :program:`poly_to_txt` tool converts :ref:`polytropic stellar
models <poly-models>` in the :ref:`POLY <poly-file-format>` format, to
a simple text-based format where the first line specifies the column
names, and subsequent lines contain the column rows.

Options
-------

.. option:: -h, --help

   Print a summary of options.

.. option:: --drop_outer

   Drop the outer-most point (which is singular when the surface
   pressure/density vanish).
