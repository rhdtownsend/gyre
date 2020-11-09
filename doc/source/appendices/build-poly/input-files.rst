.. _build-poly-input:

Input Files
===========

The :program:`build_poly` executable reads parameters from an input
file that defines a number of Fortran namelist groups, as described
below.

Polytrope Parameters
--------------------

The :nml_g:`poly` namelist group defines polytrope parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`n_r` (default :nml_v:`1`)
  Number of regions

:nml_n:`n_poly` (default :nml_v:`0`)
  Comma-sepratated list of length :nml_n:`n_r`, specifying polytropic indices for regions

:nml_n:`z_b`
  Comma-separated list of length :nml_n:`n_r`-1, specifying radial coordinates of boundaries
  between regions

:nml_n:`Delta_b`
  Comma-separated list of length :nml_n:`n_r`-1, specifying logarithmic density jumps at boundaries
  between regions

:nml_n:`Gamma_1` (default :nml_v:`5./3.`)
  First adiabatic exponent

Numerical Parameters
--------------------

The :nml_g:`num` namelist group defines numerical parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`dz` (default :nml_v:`1E-2`)
  Spacing of grid points in polytropic radial coordinate :math:`z`

:nml_n:`toler` (default :nml_v:`1E-10`)
  Relative and absolute tolerance of Lane-Emden integrations

Output Parameters
-----------------

The :nml_g:`out` namelist group defines output parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`file`
  Name of :ref:`POLY-format <poly-file-format>` file to write to
  
