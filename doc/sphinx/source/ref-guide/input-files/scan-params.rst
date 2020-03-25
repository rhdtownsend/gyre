.. _scan-params:

Frequency Scan Parameters
=========================

The :nml_g:`scan` namelist group(s) defines frequency scan parameters;
the input file can contain one or more, and points defined by each
:nml_g:`scan` namelist group are merged together. Allowable fields
are:

:nml_n:`grid_type` (default :nml_v:`'LINEAR'`)
  Distribution of frequency points; one of:

  - :nml_v:`'LINEAR'` : Uniform in frequency
  - :nml_v:`'INVERSE'` : Uniform in inverse frequency (i.e., period)

:nml_n:`grid_frame` (default :nml_v:`'INERTIAL'`)
  Frame in which :nml_n:`grid_type` applies; one of:

  - :nml_v:`'INERTIAL'` : Inertial frame
  - :nml_v:`'COROT_I'` : Co-rotating frame at inner boundary
  - :nml_v:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_n:`freq_min` (default :nml_v:`1`)
  Minimum frequency

:nml_n:`freq_max` (default :nml_v:`10`)
  Maximum frequency
  
:nml_n:`n_freq` (default :nml_v:`10`)
  Number of frequency points
  
:nml_n:`freq_min_units` (default :nml_v:`'NONE'`)
  Units of :nml_n:`freq_min`; one of:

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : linear frequency in Hz [#only_evol]_
  - :nml_v:`'UHZ'` : linear frequency in μHz [#only_evol]_
  - :nml_v:`'RAD_PER_SEC'` : angular frequency in radians per second [#only_evol]_
  - :nml_v:`'CYC_PER_DAY'` : linear frequency in cycles per day [#only_evol]_
  - :nml_v:`'ACOUSTIC_DELTA'` : Fraction of the asymptotic acoustic large frequency separation :math:`\Delta \nu`
  - :nml_v:`'GRAVITY_DELTA'` : Fraction of the asymptotic inverse gravity period separation :math:`(\Delta P)^{-1}`
  - :nml_v:`'UPPER_DELTA'` : Greater of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_v:`'LOWER_DELTA'` : Lesser of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_v:`'ACOUSTIC_CUTOFF'` : fraction of the acoustic cutoff frequency [#only_evol]_
  - :nml_v:`'GRAVITY_CUTOFF'` : fraction of the gravity cutoff frequency [#only_evol]_
  - :nml_v:`'ROSSBY_I'` : fraction of Rossby frequency at inner boundary
  - :nml_v:`'ROSSBY_O'` : fraction of Rossby frequency at outer boundary

:nml_n:`freq_max_units` (default :nml_v:`'NONE'`)
  Units of :nml_n:`freq_max`; same options as :nml_n:`freq_min_units`
  
:nml_n:`freq_min_frame` (default :nml_v:`'INERTIAL'`)
  Frame of :nml_n:`freq_min`; one of:

   - :nml_v:`'INERTIAL'` : Inertial frame
   - :nml_v:`'COROT_I'` : Co-rotating frame at inner boundary
   - :nml_v:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_n:`freq_max_frame` (default :nml_v:`'INERTIAL'`)
  Frame of :nml_n:`freq_max`; same options as :nml_n:`freq_min_frame`
  
:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_n:`model_type` is :nml_v:`'EVOL'`
