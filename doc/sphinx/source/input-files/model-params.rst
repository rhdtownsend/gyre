Stellar Model Parameters
========================

The :nml_g:`model` namelist group defines stellar model parameters;
the input file should contain exactly one. Allowable fields are:

:nml_o:`model_type`
  Type of model to use; one of:

  - :nml_l:`'HOM'` : Homogeneous compressible model
  - :nml_l:`'POLY'` : Polytropic model read from external file
  - :nml_l:`'EVOL'` : Evolutionary model read from external file

:nml_o:`file`
  Name of file, when :nml_o:`model_type` is :nml_l:`'POLY'` or :nml_l:`'EVOL'`

:nml_o:`file_format`
  Format of file, when :nml_o:`model_type` is :nml_l:`'EVOL'`; one of

  - :nml_l:`'AMDL'` : AMDL-format binary file
  - :nml_l:`'B3'` : B3-format HDF5 file
  - :nml_l:`'FAMDL'` : FAMDL-format text file
  - :nml_l:`'FGONG'` : FGONG-format text file
  - :nml_l:`'GSM'` : GSM-format HDF5 file
  - :nml_l:`'LOSC'` : LOSC-format text file
  - :nml_l:`'MESA'` : MESA GYRE-format text file
  - :nml_l:`'OSC'` : OSC-format text file

:nml_o:`data_format` (default :nml_l:`''`, indicates auto-select)
  Fortran format specifier for data read from OSC-, FGONG- and FAMDL-format files
  
:nml_o:`deriv_type` (default :nml_l:`'MONO'`)
  Cubic interpolation derivatives type, when :nml_o:`model_type` is :nml_l:`'POLY'` or :nml_l:`'EVOL'`; one of

  - :nml_l:`'NATURAL'` : Natural (spline) derivatives
  - :nml_l:`'FINDIFF'` : Finite-difference derivatives
  - :nml_l:`'MONO'` : Monotonized derivatives (default)

:nml_o:`Gamma_1` (default :nml_l:`5/3`)
  First adiabatic exponent, when :nml_o:`model_type` is :nml_l:`'HOM'`
   
:nml_o:`grid_type` (default :nml_l:`'UNI'`)
  Model grid type, when :nml_o:`model_type` is :nml_l:`'HOM'`; one of

  - :nml_l:`'UNI'` : Uniform spacing
  - :nml_l:`'GEO'` : Geometric spacing
  - :nml_l:`'LOG'` : Logarithmic spacing

:nml_o:`n` (default :nml_l:`10`)
  Number of points in model grid, when :nml_o:`model_type` is :nml_l:`'HOM'`
       
:nml_o:`s` (default :nml_l:`1`)
  Skewness parameter for model grid, when
  :nml_o:`model_type` is :nml_l:`'HOM'` and :nml_o:`grid_type` is
  :nml_l:`'GEO'` or :nml_l:`'LOG'`

:nml_o:`x_i` (default :nml_l:`0`)
  Inner boundary coordinate of model grid, when :nml_o:`model_type` is :nml_l:`'HOM'`
    
:nml_o:`x_o` (default :nml_l:`1`)
  Outer boundary coordinate of model grid, when :nml_o:`model_type` is :nml_l:`'HOM'`

:nml_o:`uniform_rot` (default :nml_l:`.FALSE.`)
  Flag to force uniform rotation

:nml_o:`Omega_rot` (default :nml_l:`.FALSE.`)
  Rotation angular velocity, when :nml_o:`uniform_rot` is :nml_l:`.TRUE.`

:nml_o:`Omega_units` (default :nml_l:`'NONE'`)
  Units of :nml_o:`Omega_rot`; one of

  - :nml_l:`'NONE'` : Dimensionless angular frequency
  - :nml_l:`'HZ'` : Linear frequency in Hz [#only_evol]_
  - :nml_l:`'UHZ'` : Linear frequency in μHz [#only_evol]_
  - :nml_l:`'RAD_PER_SEC'` : Angular frequency in radians per second [#only_evol]_
  - :nml_l:`'CYC_PER_DAY'` : Linear frequency in cycles per day [#only_evol]_
  - :nml_l:`'CRITICAL'` : Fraction of the Roche critical rate [#only_evol]_

:nml_o:`dx_snap` (default :nml_l:`0`)
  Threshold for snapping model points together, when
  :nml_o:`model_type` is :nml_l:`'EVOL'`. If a pair of points are
  separated by less than :nml_l:`dx_snap`, they are snapped together.

:nml_o:`add_center` (default :nml_l:`.TRUE.`)
  Flag to add a center point to the model, when :nml_o:`model_type` is
  :nml_l:`'EVOL'` or :nml_l:`'POLY'`. If a point does not already
  exist at the origin, then one is added

:nml_o:`repair_As` (default :nml_l:`.FALSE.`)
  Flag to repair inaccuracies in the dimensionless Brunt-Väisälä
  frequency at density discontinuities

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_o:`model_type` is :nml_l:`'EVOL'`
