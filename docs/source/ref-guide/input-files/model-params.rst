.. _model-params:

Stellar Model Parameters
========================

The :nml_g:`model` namelist group defines stellar model parameters, as
follows:

:nml_n:`model_type`
  Type of model to use; one of:

  - :nml_v:`'HOM'` : :ref:`Homogeneous compressible model <hom-models>`
  - :nml_v:`'POLY'` : :ref:`Polytropic model <poly-models>` read from external file
  - :nml_v:`'ANAPOLY'` : :ref:`Analytic polytropic model <anapoly-models>`
  - :nml_v:`'EVOL'` : :ref:`Evolutionary model <evol-models>` read from external file

:nml_n:`file`
  Name of file (when :nml_n:`model_type`\ = \ :nml_v:`'POLY'`\ \|\ :nml_v:`'EVOL'`)

:nml_n:`file_format`
  Format of file (when :nml_n:`model_type`\ = \ :nml_v:`'EVOL'`); one of

  - :nml_v:`'AMDL'` : AMDL-format binary file
  - :nml_v:`'B3'` : B3-format HDF5 file
  - :nml_v:`'FAMDL'` : FAMDL-format text file
  - :nml_v:`'FGONG'` : FGONG-format text file
  - :nml_v:`'GSM'` : :ref:`GSM-format <gsm-file-format>` HDF5 file
  - :nml_v:`'LOSC'` : LOSC-format text file
  - :nml_v:`'MESA'` : :ref:`MESA/GYRE-format <mesa-file-format>` text file
  - :nml_v:`'OSC'` : OSC-format text file
  - :nml_v:`'POLY'` : :ref:`POLY-format <poly-file-format>` HDF5 file
  - :nml_v:`'WDEC'` : WDEC-format text file

  Using :nml_n:`model_type`\ =\ :nml_v:`'EVOL'` and
  :nml_n:`file_format`\ =\ :nml_v:`'POLY'` gives similar (but not
  identical) results to choosing :nml_n:`model_type`\ =\
  :nml_v:`'POLY'`; the two approaches differ in how the model
  structure is interpolated.

:nml_n:`data_format` (default :nml_v:`''`, indicates auto-select)
  Fortran format specifier for data read from OSC-, FGONG- and FAMDL-format files
  
:nml_n:`deriv_type` (default :nml_v:`'MONO'`)
  Cubic interpolation derivatives type (when :nml_n:`model_type`\ =\ :nml_v:`'EVOL'`); one of

  - :nml_v:`'SPLINE'` : Spline (non-local) derivatives
  - :nml_v:`'FINDIFF'` : Finite-difference derivatives
  - :nml_v:`'MONO'` : Monotonized derivatives (default)

:nml_n:`Gamma_1` (default :nml_v:`5/3`)
  First adiabatic exponent (when :nml_n:`model_type`\ =\ :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'`)

:nml_n:`n_poly` (default :nml_n:`0`)
  Polytropic index (when :nml_n:`model_type`\ =\ :nml_v:`'ANAPOLY'`)

:nml_n:`theta_s` (default :nml_n:`0`)
  Surface value of polytropic dependent variable (when :nml_n:`model_type`\ =\ :nml_v:`'ANAPOLY'`)

:nml_n:`grid_type` (default :nml_v:`'UNI'`)
  Model grid type (when :nml_n:`model_type`\ =\ :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'`); one of

  - :nml_v:`'UNI'` : Uniform spacing
  - :nml_v:`'GEO'` : Geometric spacing
  - :nml_v:`'LOG'` : Logarithmic spacing

:nml_n:`n` (default :nml_v:`10`)
  Number of points in model grid (when :nml_n:`model_type`\ =\ :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'`)
       
:nml_n:`s` (default :nml_v:`1`)
  Skewness parameter for model grid (when :nml_n:`model_type`\ =\
  :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'` and :nml_n:`grid_type`\ =\ :nml_v:`'GEO'`\ \|\
  :nml_v:`'LOG'`)

:nml_n:`x_i` (default :nml_v:`0`)
  Inner boundary coordinate of model grid (when :nml_n:`model_type`\ =\ :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'`)
    
:nml_n:`x_o` (default :nml_v:`1`)
  Outer boundary coordinate of model grid (when :nml_n:`model_type`\ =\ :nml_v:`'HOM'`\ \|\ :nml_v:`'ANAPOLY'`)

:nml_n:`dx_snap` (default :nml_v:`0`)
  Threshold for snapping model points together, when
  :nml_n:`model_type` is :nml_v:`'EVOL'`. If a pair of points are
  separated by less than :nml_v:`dx_snap`, they are snapped together.

:nml_n:`add_center` (default :nml_v:`.TRUE.`)
  Flag to add a center point to the model (when :nml_n:`model_type`\ =\
  :nml_v:`'EVOL'`\ \|\ :nml_v:`'POLY'`). If a point does not already
  exist at the origin, then one is added

:nml_n:`repair_As` (default :nml_v:`.FALSE.`)
  Flag to repair inaccuracies in the dimensionless Brunt-Väisälä
  frequency at density discontinuities

:nml_n:`repair_derivs` (default :nml_v:`.FALSE.`)
  Flag to repair derivatives of :math:`c_1` and :math:`U` structure coefficients
  (when :nml_n:`model_type`\ =\ :nml_v:`'EVOL'`), in accordance with equations
  (20) and (21) of :ads_citet:`takata:2006a`, respectively
  
