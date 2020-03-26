.. _rot-params:

Rotation Parameters
===================

The :nml_g:`rot` namelist group defines rotational parameters; the
input file can contain one or more, but only the last (matching) one
is used.  Allowable parameters are:

:nml_n:`Omega_rot` (default :nml_v:`0`)
  Rotation angular frequency

:nml_n:`coriolis_method` (default :nml_v:`'NULL'`)
  method used to treat the Coriolis force; one of:

  - :nml_v:`'NULL'` : Neglect the Coriolis force
  - :nml_v:`'TAR'` : Use the traditional approximation of rotation

:nml_n:`Omega_rot_source` (default :nml_v:`'MODEL'`)

:nml_n:`Omega_rot_units` (default :nml_v:`'NULL'`)

:nml_n:`complex_lambda` (default :nml_v:`.FALSE.`)

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
