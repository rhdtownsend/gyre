.. _rot-params:

Rotation Parameters
===================

The :nml_g:`rot` namelist group defines rotational parameters, as
follows:

:nml_n:`coriolis_method` (default :nml_v:`'NULL'`)
  Method used to treat the Coriolis force; one of:

  - :nml_v:`'NULL'` : Neglect the Coriolis force
  - :nml_v:`'TAR'` : Use the traditional approximation of rotation

:nml_n:`Omega_rot_source` (default :nml_v:`'MODEL'`)
  Source for rotational angular frequency :math:`\Orot`; one of:

  - :nml_v:`'MODEL'` : Differential rotation, with a spatially varying :math:`\Orot`
    obtained from the stellar model
  - :nml_v:`'UNIFORM'` : Uniform rotation, with a spatially constant :math:`\Orot` set
    by the :nml_n:`Omega_rot` and :nml_n:`Omega_rot_units` parameters

:nml_n:`Omega_rot` (default :nml_v:`0`)
  Rotation angular frequency (when :nml_n:`Omega_rot_source`\ =\ :nml_v:`'UNIFORM'`)

:nml_n:`Omega_rot_units` (default :nml_v:`'NONE'`)
  Units of :nml_n:`Omega_rot` (when :nml_n:`Omega_rot_source`\ =\ :nml_v:`'UNIFORM'`); one of:

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : Linear frequency in Hz\ [#only-D]_
  - :nml_v:`'UHZ'` : Linear frequency in :math:`\mu`\ Hz\ [#only-D]_
  - :nml_v:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only-D]_
  - :nml_v:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only-D]_
  - :nml_v:`'CRITICAL'` : Fraction of the Roche critical rate\ [#only-D]_

:nml_n:`rossby` (default :nml_v:`.FALSE.`)
  Flag to use Rossby solution family in TAR (when :nml_n:`coriolis_method`\ =\ :nml_v:`'TAR'`)

:nml_n:`complex_lambda` (default :nml_v:`.FALSE.`)
  Flag to use complex arithmetic when evaluating the TAR angular eigenvalue :math:`\lambda` (when :nml_n:`coriolis_method`\ =\ :nml_v:`'TAR'`)

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
  Comma-separated list of :nml_g:`mode` tags to match

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`

