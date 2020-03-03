.. _model-params:

Stellar Model Parameters
========================

The :nml_g:`model` namelist group defines stellar model parameters;
the input file should contain exactly one. Allowable fields are:

:nml_n:`model_type`
  Type of model to use; one of:

  - :nml_v:`'HOM'` : Homogeneous compressible model
  - :nml_v:`'POLY'` : Polytropic model read from external file
  - :nml_v:`'EVOL'` : Evolutionary model read from external file

:nml_n:`file`
  Name of file, when :nml_n:`model_type` is :nml_v:`'POLY'` or :nml_v:`'EVOL'`

:nml_n:`file_format`
  Format of file, when :nml_n:`model_type` is :nml_v:`'EVOL'`; one of

  - :nml_v:`'AMDL'` : AMDL-format binary file
  - :nml_v:`'B3'` : B3-format HDF5 file
  - :nml_v:`'FAMDL'` : FAMDL-format text file
  - :nml_v:`'FGONG'` : FGONG-format text file
  - :nml_v:`'GSM'` : GSM-format HDF5 file
  - :nml_v:`'LOSC'` : LOSC-format text file
  - :nml_v:`'MESA'` : MESA GYRE-format text file
  - :nml_v:`'OSC'` : OSC-format text file

:nml_n:`data_format` (default :nml_v:`''`, indicates auto-select)
  Fortran format specifier for data read from OSC-, FGONG- and FAMDL-format files
  
:nml_n:`deriv_type` (default :nml_v:`'MONO'`)
  Cubic interpolation derivatives type, when :nml_n:`model_type` is :nml_v:`'POLY'` or :nml_v:`'EVOL'`; one of

  - :nml_v:`'NATURAL'` : Natural (spline) derivatives
  - :nml_v:`'FINDIFF'` : Finite-difference derivatives
  - :nml_v:`'MONO'` : Monotonized derivatives (default)

:nml_n:`Gamma_1` (default :nml_v:`5/3`)
  First adiabatic exponent, when :nml_n:`model_type` is :nml_v:`'HOM'`
   
:nml_n:`grid_type` (default :nml_v:`'UNI'`)
  Model grid type, when :nml_n:`model_type` is :nml_v:`'HOM'`; one of

  - :nml_v:`'UNI'` : Uniform spacing
  - :nml_v:`'GEO'` : Geometric spacing
  - :nml_v:`'LOG'` : Logarithmic spacing

:nml_n:`n` (default :nml_v:`10`)
  Number of points in model grid, when :nml_n:`model_type` is :nml_v:`'HOM'`
       
:nml_n:`s` (default :nml_v:`1`)
  Skewness parameter for model grid, when
  :nml_n:`model_type` is :nml_v:`'HOM'` and :nml_n:`grid_type` is
  :nml_v:`'GEO'` or :nml_v:`'LOG'`

:nml_n:`x_i` (default :nml_v:`0`)
  Inner boundary coordinate of model grid, when :nml_n:`model_type` is :nml_v:`'HOM'`
    
:nml_n:`x_o` (default :nml_v:`1`)
  Outer boundary coordinate of model grid, when :nml_n:`model_type` is :nml_v:`'HOM'`

:nml_n:`uniform_rot` (default :nml_v:`.FALSE.`)
  Flag to force uniform rotation

:nml_n:`Omega_rot` (default :nml_v:`.FALSE.`)
  Rotation angular velocity, when :nml_n:`uniform_rot` is :nml_v:`.TRUE.`

:nml_n:`Omega_units` (default :nml_v:`'NONE'`)
  Units of :nml_n:`Omega_rot`; one of

  - :nml_v:`'NONE'` : Dimensionless angular frequency
  - :nml_v:`'HZ'` : Linear frequency in Hz [#only_evol]_
  - :nml_v:`'UHZ'` : Linear frequency in μHz [#only_evol]_
  - :nml_v:`'RAD_PER_SEC'` : Angular frequency in radians per second [#only_evol]_
  - :nml_v:`'CYC_PER_DAY'` : Linear frequency in cycles per day [#only_evol]_
  - :nml_v:`'CRITICAL'` : Fraction of the Roche critical rate [#only_evol]_

:nml_n:`dx_snap` (default :nml_v:`0`)
  Threshold for snapping model points together, when
  :nml_n:`model_type` is :nml_v:`'EVOL'`. If a pair of points are
  separated by less than :nml_v:`dx_snap`, they are snapped together.

:nml_n:`add_center` (default :nml_v:`.TRUE.`)
  Flag to add a center point to the model, when :nml_n:`model_type` is
  :nml_v:`'EVOL'` or :nml_v:`'POLY'`. If a point does not already
  exist at the origin, then one is added

:nml_n:`repair_As` (default :nml_v:`.FALSE.`)
  Flag to repair inaccuracies in the dimensionless Brunt-Väisälä
  frequency at density discontinuities

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_n:`model_type` is :nml_v:`'EVOL'`
