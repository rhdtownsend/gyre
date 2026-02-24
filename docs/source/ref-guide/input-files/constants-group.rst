.. _constants-group:

.. nml:group:: constants

Constants Namelist Group
========================

The :nml:group:`constants` namelist group defines various physical constants. The following options are available:

.. nml:option:: l
   :type: real

   Gravitational constant :math:`G`

.. nml:option:: C_LIGHT
   :type: real

   Speed of light *in vacuo* :math:`c`

.. nml:option:: A_RADIATION
   :type: real

   Radiation constant :math:`a`

.. nml:option:: M_SUN
   :type: real

   Solar mass :math:`\Msun`

.. nml:option:: R_SUN
   :type: real

   Solar radius :math:`\Rsun`

.. nml:option:: L_SUN
   :type: real

   Solar luminosity :math:`\Lsun`

.. nml:option:: GYRE_DIR
   :type: string

   Top-level GYRE directory; overrides the :envvar:`GYRE_DIR`
   environment variable

All of these constants are in cgs units (where applicable), and the
default values are defined in
:file:`{$GYRE_DIR}/src/common/gyre_constants.fpp`.
