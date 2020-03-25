.. _constants:

Constants
=========

The :nml_g:`constants` namelist group defines various physical
constants; the input file should contain exactly one. Allowable fields
are:

:nml_n:`G_GRAVITY`
    Gravitational constant :math:`G`

:nml_n:`C_LIGHT`
    Speed of light *in vacuo* :math:`c`

:nml_n:`A_RADIATION`
    Radiation constant :math:`a`

:nml_n:`M_SUN`
    Solar mass :math:`M_{\odot}`

:nml_n:`R_SUN`
    Solar radius :math:`R_{\odot}`

:nml_n:`L_SUN`
    Solar luminosity :math:`L_{\odot}`

:nml_n:`GYRE_DIR`
    Top-level GYRE directory; overrides the ``GYRE_DIR``
    environment variable

All of these constants are in cgs units (where applicable), and the
default values are defined in :repo:`gyre_constants.fpp
<src/common/gyre_constants.fpp>`.
