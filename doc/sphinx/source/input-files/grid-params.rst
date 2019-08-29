Grid Parameters
===============

The :nml_g:`grid` namelist group defines the parameters used to set up
the calculation grid; the input file can contain one or more, but only
the last (matching) one is used. Allowable fields are:

:nml_o:`x_i` (default based on model grid)
  Inner boundary coordinate of calculation grid

:nml_o:`x_o` (default based on model grid)
  Outer boundary coordinate of calculation grid

:nml_o:`n_inner` (default :nml_l:`0`)
  Number of extra point to insert at center

:nml_o:`alpha_osc` (default :nml_l:`0`)
  Oscillatory resampling parameter :math:`\alpha_{\rm osc}`

:nml_o:`alpha_exp` (default :nml_l:`0`)
  Exponential resampling parameter :math:`\alpha_{\rm exp}`

:nml_o:`alpha_thm` (default :nml_l:`0`)
  Thermal resampling parameter :math:`\alpha_{\rm thm}`

:nml_o:`alpha_str` (default :nml_l:`0`)
  Structural resampling parameter :math:`\alpha_{\rm str}`

:nml_o:`dx_min` (default :nml_l:`SQRT(EPSILON(1._WP))`)
  Minimum spacing of grid points
  
:nml_o:`n_iter_max` (default :nml_l:`8`)
  Maximum number of times to iteratively resample

:nml_o:`tag_list` (default :nml_l:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

See `Understanding Grids <Understanding Grids>`__ for further details,
in particular a discussion of how the :nml_o:`alpha_*` resampling
parameters work.
