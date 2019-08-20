Frequency Scan Parameters
=========================

The :nml_g:`scan` namelist group(s) defines frequency scan parameters;
the input file can contain one or more, and points defined by each
:nml_g:`scan` namelist group are merged together. Allowable fields
are:

:nml_o:`grid_type` (default :nml_l:`'LINEAR'`)
  Distribution of frequency points; one of:

  - :nml_l:`'LINEAR'` : Uniform in frequency
  - :nml_l:`'INVERSE'` : Uniform in inverse frequency (i.e., period)

:nml_o:`grid_frame` (default :nml_l:`'INERTIAL'`)
  Frame in which :nml_o:`grid_type` applies; one of:

  - :nml_l:`'INERTIAL'` : Inertial frame
  - :nml_l:`'COROT_I'` : Co-rotating frame at inner boundary
  - :nml_l:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_o:`freq_min` (default :nml_l:`1`)
  Minimum frequency

:nml_o:`freq_max` (default :nml_l:`10`)
  Maximum frequency
  
:nml_o:`n_freq` (default :nml_l:`10`)
  Number of frequency points
  
:nml_o:`freq_min_units` (default :nml_l:`'NONE'`)
  Units of :nml_o:`freq_min`; one of:

  - :nml_l:`'NONE'` : Dimensionless angular frequency
  - :nml_l:`'HZ'` : linear frequency in Hz [#only_evol]_
  - :nml_l:`'UHZ'` : linear frequency in μHz [#only_evol]_
  - :nml_l:`'RAD_PER_SEC'` : angular frequency in radians per second [#only_evol]_
  - :nml_l:`'CYC_PER_DAY'` : linear frequency in cycles per day [#only_evol]_
  - :nml_l:`'ACOUSTIC_DELTA'` : Fraction of the acoustic large frequency separation :math:`\Delta \nu`
  - :nml_l:`'GRAVITY_DELTA'` : Fraction of the inverse gravity period separation :math:`(\Delta P)^{-1}`
  - :nml_l:`'UPPER_DELTA'` : Greater of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_l:`'LOWER_DELTA'` : Lesser of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_l:`'ACOUSTIC_CUTOFF'` : fraction of the acoustic cutoff frequency [#only_evol]_
  - :nml_l:`'GRAVITY_CUTOFF'` : fraction of the gravity cutoff frequency [#only_evol]_
  - :nml_l:`'ROSSBY_I'` : fraction of Rossby frequency at inner boundary
  - :nml_l:`'ROSSBY_O'` : fraction of Rossby frequency at outer boundary

:nml_o:`freq_max_units` (default :nml_l:`'NONE'`)
  Units of :nml_o:`freq_max`; same options as :nml_o:`freq_min_units`
  
:nml_o:`freq_min_frame` (default :nml_l:`'INERTIAL'`)
  Frame of :nml_o:`freq_min`; one of:

   - :nml_l:`'INERTIAL'` : Inertial frame
   - :nml_l:`'COROT_I'` : Co-rotating frame at inner boundary
   - :nml_l:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_o:`freq_max_frame` (default :nml_l:`'INERTIAL'`)
  Frame of :nml_o:`freq_max`; same options as :nml_o:`freq_min_frame`
  
:nml_o:`tag_list` (default :nml_l:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_o:`model_type` is :nml_l:`'EVOL'`
