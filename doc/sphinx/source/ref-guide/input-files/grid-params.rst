.. _grid-params:

Grid Parameters
===============

The :nml_g:`grid` namelist group defines the parameters used to set up
the calculation grid; the input file can contain one or more, but only
the last (matching) one is used. Allowable fields are:

:nml_n:`x_i` (default based on model grid)
  Inner boundary coordinate of calculation grid

:nml_n:`x_o` (default based on model grid)
  Outer boundary coordinate of calculation grid

:nml_n:`n_inner` (default :nml_v:`0`)
  Number of extra point to insert at center

:nml_n:`alpha_osc` (default :nml_v:`0`)
  Oscillatory resampling parameter :math:`\alpha_{\rm osc}`

:nml_n:`alpha_exp` (default :nml_v:`0`)
  Exponential resampling parameter :math:`\alpha_{\rm exp}`

:nml_n:`alpha_thm` (default :nml_v:`0`)
  Thermal resampling parameter :math:`\alpha_{\rm thm}`

:nml_n:`alpha_str` (default :nml_v:`0`)
  Structural resampling parameter :math:`\alpha_{\rm str}`

:nml_n:`dx_min` (default :nml_v:`SQRT(EPSILON(1._WP))`)
  Minimum spacing of grid points
  
:nml_n:`n_iter_max` (default :nml_v:`8`)
  Maximum number of times to iteratively resample

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

See `Understanding Grids <Understanding Grids>`__ for further details,
in particular a discussion of how the :nml_n:`alpha_*` resampling
parameters work.
