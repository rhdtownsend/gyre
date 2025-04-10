.. _support-tools-poly-to-txt:

poly_to_txt
===========

.. program:: poly_to_txt

The :program:`poly_to_txt` tool converts :ref:`polytropic stellar
models <poly-models>` in the :ref:`POLY <poly-file-format>` format, to
a simple text-based format.

Syntax
------

:program:`poly_to_txt` :option:`in_filename` :option:`out_filename` :option:`drop_outer`

Parameters
----------

.. option:: in_filename

   File name of input POLY file

.. option:: out_filename

   File name of output text file

.. option:: drop_outer

   Flag to drop outer-most point (which is singular when the surface pressure/density vanish)

Output
------

The output is written to the :option:`out_filename` file. This file is
a text (ASCII) table. The first line specify the column names, and
subsequent lines give the column rows.
