.. _rot-params:

Rotation Parameters
===================

The :nml_g:`rot` namelist group defines rotational parameters; the
input file can contain one or more, but only the last (matching) one
is used.  Allowable parameters are:

:nml_n:`coriolis_method` (default :nml_v:`'NULL'`)
  Method used to treat the Coriolis force; one of:

  - :nml_v:`'NULL'` : Neglect the Coriolis force
  - :nml_v:`'TAR'` : Use the traditional approximation of rotation

:nml_n:`Omega_rot_source` (default :nml_v:`'MODEL'`)
  Sourcce for rotational angular velocity; one of:

  - :nml_v:`'MODEL'` : Use data from the stellar model
  - :nml_v:`'UNIFORM'` : Uniform rotation, with angular velocity set by :nml_n:`Omega_rot` parameter

:nml_n:`Omega_rot` (default :nml_v:`0`)
  Rotation angular frequency (:nml_n:`Omega_rot_source`\ ==\ :nml_v:`'UNIFORM'`)

:nml_n:`Omega_rot_units` (default :nml_v:`'NULL'`)
  Units of :nml_n:`Omega_rot` (:nml_n:`Omega_rot_source`\ ==\ :nml_v:`'UNIFORM'`); one of:

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : Linear frequency in Hz\ [#only_evol]_
  - :nml_v:`'UHZ'` : Linear frequency in μHz\ [#only_evol]_
  - :nml_v:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only_evol]_
  - :nml_v:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only_evol]_
  - :nml_v:`'CRITICAL'` : Fraction of the Roche critical rate\ [#only_evol]_

:nml_n:`rossby` (default :nml_v:`.FALSE.`)
  Flag to use Rossby solution family in TAR

:nml_n:`complex_lambda` (default :nml_v:`.FALSE.`)
  Flag to use complex arithmetic when evaluating the TAR angular eigenvalue :math:`\lambda`

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
  Comma-separated list of :nml_g:`mode` tags to match

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_n:`model_type`\ ==\ :nml_v:`'EVOL'`
  
