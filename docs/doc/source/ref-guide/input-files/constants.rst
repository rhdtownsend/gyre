.. _constants:

Constants
=========

The :nml_g:`constants` namelist group defines various physical
constants, as follows:

:nml_n:`G_GRAVITY`
    Gravitational constant :math:`G`

:nml_n:`C_LIGHT`
    Speed of light *in vacuo* :math:`c`

:nml_n:`A_RADIATION`
    Radiation constant :math:`a`

:nml_n:`M_SUN`
    Solar mass :math:`\Msun`

:nml_n:`R_SUN`
    Solar radius :math:`\Rsun`

:nml_n:`L_SUN`
    Solar luminosity :math:`\Lsun`

:nml_n:`GYRE_DIR`
    Top-level GYRE directory; overrides the :envvar:`GYRE_DIR`
    environment variable

All of these constants are in cgs units (where applicable), and the
default values are defined in
:file:`{$GYRE_DIR}/src/common/gyre_constants.fpp`.
