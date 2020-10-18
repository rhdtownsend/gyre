.. _comp-ptrope-input:

Input Files
===========

The :program:`build_poly` executable reads parameters from an input
file that defines a number of Fortran namelist groups, as described
below.

Polytrope Parameters
--------------------

The :nml_g:`poly` namelist group defines polytrope parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`n_d` (default :nml_v:`0`)
  Number of points at which the density and polytropic index can jump;
  establishes the length of the :nml_n:`n_poly`, :nml_n:`xi_d` and
  :nml_n:`Delta_d` arrays

:nml_n:`n_poly` (default :nml_v:`0`)
  Comma-sepratated list of length :nml_n:`n_d`+1, specifying polytropic indices for regions

:nml_n:`xi_d`
  Comma-separated list of length :nml_n:`n_d`, specifying coordinates of boundaries
  between regions

:nml_n:`Delta_d`
  Comma-separated list of length :nml_n:`n_d`, specifying logarithmic density jumps at boundaries
  between regions

:nml_n:`Gamma_1` (default :nml_v:`1.666666...`)
  First adiabatic exponent

Numerical Parameters
--------------------

The :nml_g:`num` namelist group defines numerical parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`dxi` (default :nml_v:``)
  Spacing of grid points in polytropic radial coordinate :math:`\xi`

:nml_n:`toler` (default :nml_v:``)
  Relative tolerance of Lane-Emden integrations

Output Parameters
-----------------

The :nml_g:`out` namelist group defines output parameters; the
input file can contain only one. Allowable parameters are:

:nml_n:`file`
  Name of file to write `POLY`-format file to
  
  

