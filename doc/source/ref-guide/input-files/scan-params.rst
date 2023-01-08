.. _scan-params:

Frequency Scan Parameters
=========================

The :nml_g:`scan` namelist group defines frequency grid parameters, as
follows:

:nml_n:`grid_type` (default :nml_v:`'LINEAR'`)
  Distribution of frequency points; one of:

  - :nml_v:`'LINEAR'` : Uniform in frequency
  - :nml_v:`'INVERSE'` : Uniform in inverse frequency
  - :nml_v:`'FILE'` : Read from file

:nml_n:`grid_frame` (default :nml_v:`'INERTIAL'`)
  Reference frame in which :nml_n:`grid_type` applies; one of:

  - :nml_v:`'INERTIAL'` : Inertial frame
  - :nml_v:`'COROT_I'` : Co-rotating frame at inner boundary
  - :nml_v:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_n:`freq_min` (default :nml_v:`1`)
  Minimum frequency, when :nml_n:`grid_type` is :nml_v:`'LINEAR'` or :nml_v:`'INVERSE'`

:nml_n:`freq_max` (default :nml_v:`10`)
  Maximum frequency, when :nml_n:`grid_type` is :nml_v:`'LINEAR'` or :nml_v:`'INVERSE'`
  
:nml_n:`n_freq` (default :nml_v:`10`)
  Number of frequency points, when :nml_n:`grid_type` is :nml_v:`'LINEAR'` or :nml_v:`'INVERSE'`

:nml_n:`freq_units` (default :nml_v:`NONE`)
  Units of :nml_n:`freq_min` and :nml_n:`freq_max`, when
  :nml_n:`grid_type` is :nml_v:`'LINEAR'` or :nml_v:`'INVERSE'`; units
  of read frequencies when :nml_n:`grid_type` is :nml_v:`'FILE'`

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : linear frequency in Hz\ [#only-D]_
  - :nml_v:`'UHZ'` : linear frequency in μHz\ [#only-D]_
  - :nml_v:`'RAD_PER_SEC'` : angular frequency in radians per second\ [#only-D]_
  - :nml_v:`'CYC_PER_DAY'` : linear frequency in cycles per day\ [#only-D]_
  - :nml_v:`'ACOUSTIC_DELTA'` : Fraction of the asymptotic acoustic large frequency separation :math:`\Delta \nu`
  - :nml_v:`'GRAVITY_DELTA'` : Fraction of the asymptotic inverse gravity period separation :math:`(\Delta P)^{-1}`
  - :nml_v:`'UPPER_DELTA'` : Greater of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_v:`'LOWER_DELTA'` : Lesser of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_v:`'ACOUSTIC_CUTOFF'` : fraction of the acoustic cutoff frequency\ [#only-D]_
  - :nml_v:`'GRAVITY_CUTOFF'` : fraction of the gravity cutoff frequency\ [#only-D]_
  - :nml_v:`'ROSSBY_I'` : fraction of Rossby frequency (see equation :eq:`ross-freq`) at inner boundary
  - :nml_v:`'ROSSBY_O'` : fraction of Rossby frequency (see equation :eq:`ross-freq`) at outer boundary

:nml_n:`freq_min_units` (default :nml_v:`''`)
  Units of :nml_n:`freq_min`; same options as :nml_n:`freq_units` and overrides it if set

:nml_n:`freq_max_units` (default :nml_v:`''`)
  Units of :nml_n:`freq_max`; same options as :nml_n:`freq_units` and overrides it if set

:nml_n:`freq_frame` (default :nml_v:`'INERTIAL'`)
  Reference frame in which :nml_n:`freq_min` and :nml_n:`freq_max` are defined, when :nml_n:`grid_type`
  is :nml_v:`'LINEAR'` or :nml_v:`'INVERSE'`; one of:

   - :nml_v:`'INERTIAL'` : Inertial frame
   - :nml_v:`'COROT_I'` : Co-rotating frame at inner boundary
   - :nml_v:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_n:`file`
  File to read frequencies from, when :nml_n:`grid_type` is :nml_v:`'FILE'`

:nml_n:`axis` (default :nml_v:`'REAL`')
  Axis that scan applies to; one of

  - :nml_v:`'REAL'` : Real axis
  - :nml_v:`'IMAG'` : Imaginary axis

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

An input file can contain one or more :nml_g:`scan` namelist group;
the points defined by each (tag-matching) group are merged together to
build the frequency grid. See the :ref:`freq-grids` section for
further details.

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
