.. _orbit-params:

Orbital Parameters
==================

The :nml_g:`orbit` namelist group defines orbital
parameters, as follows:

:nml_n:`Omega_orb` (default :nml_v:`1`)
  Orbital angular frequency

:nml_n:`Omega_orb_units` (default :nml_v:`'NULL'`)
  Units of :nml_n:`Omega_orb`; one of:

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : Linear frequency in Hz\ [#only-D]_
  - :nml_v:`'UHZ'` : Linear frequency in μHz\ [#only-D]_
  - :nml_v:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only-D]_
  - :nml_v:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only-D]_
  - :nml_v:`'CRITICAL'` : Fraction of the Roche critical rate\ [#only-D]_

:nml_n:`q` (default :nml_v:`1`)
   Ratio of secondary mass to primary mass

:nml_n:`e` (default :nml_v:`0`)
   Orbital eccentricity

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
  Comma-separated list of :nml_g:`tide` tags to match

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
