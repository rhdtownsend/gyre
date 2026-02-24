.. _mode-group:

.. nml:group:: mode

Mode Namelist Group
===================

The :nml:group:`model` namelist group defines mode indices. The
following options are available:

.. nml:option:: l
   :type: integer
   :default: 0

   Harmonic degree :math:`\ell`

.. nml:option:: m
   :type: integer
   :default: 0

   Azimuthal order :math:`m`

.. nml:option:: n_pg_min
   :type: integer
   :default: -HUGE

   Filter for minimum radial order (applies only to adiabatic calculations)

.. nml:option:: n_pg_max
   :type: integer
   :default: +HUGE

   Filter for maximum radial order (applies only to adiabatic calculations)

.. nml:option:: tag
   :type: string

   Tag for controlling selection of other namelist groups (see
   :ref:`working-with-tags`)
