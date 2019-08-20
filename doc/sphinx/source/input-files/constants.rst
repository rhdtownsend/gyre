Constants
=========

The :nml_g:`constants` namelist group defines various physical
constants; the input file should contain exactly one. Allowable fields
are:

:nml_o:`G_GRAVITY`
    Gravitational constant

:nml_o:`C_LIGHT`
    Speed of light *in vacuo*

:nml_o:`A_RADIATION`
    Radiation constant

:nml_o:`M_SUN`
    Solar mass

:nml_o:`R_SUN`
    Solar radius

:nml_o:`L_SUN`
    Solar luminosity

:nml_o:`GYRE_DIR`
    Top-level GYRE directory (overrides the ``GYRE_DIR``
    environment variable)

All of these constants are in cgs units (where applicable), and the
default values are defined in :repo:`gyre_constants.fpp <src/common/gyre_constants.fpp>`.
