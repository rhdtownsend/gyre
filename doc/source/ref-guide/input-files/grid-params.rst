.. _grid-params:

Grid Parameters
===============

The :nml_g:`grid` namelist group defines the parameters used to set up
the spatial grid, as follows:

:nml_n:`grid_source` (default :nml_v:`MODEL`)
  Source for skeleton grid; one of:

  - :nml_v:`'MODEL'` : Obtained from the stellar model
  - :nml_v:`'DETAIL_FILE` : Read from a detail file, which must contain
    the :nml_v:`x` item
  - :nml_v:`'TEXT_FILE` : Read from a text file, with one :math:`x` per line

:nml_n:`x_i` (default based on model grid)
  Inner boundary coordinate of calculation grid

:nml_n:`x_o` (default based on model grid)
  Outer boundary coordinate of calculation grid

:nml_n:`w_osc` (default :nml_v:`0`)
  Oscillatory weighting parameter :math:`w_{\rm osc}`

:nml_n:`w_exp` (default :nml_v:`0`)
  Exponential weighting parameter :math:`w_{\rm exp}`

:nml_n:`w_ctr` (default :nml_v:`0`)
  Center weighting parameter :math:`w_{\rm ctr}`

:nml_n:`w_thm` (default :nml_v:`0`)
  Thermal weighting parameter :math:`w_{\rm thm}`

:nml_n:`w_str` (default :nml_v:`0`)
  Structural weighting parameter :math:`w_{\rm str}`

:nml_n:`dx_min` (default :nml_v:`SQRT(EPSILON(1._WP))`)
  Minimum spacing of grid points
  
:nml_n:`dx_max` (default :nml_v:`HUGE(0._WP)`)
  Maximum spacing of grid points
  
:nml_n:`n_iter_max` (default :nml_v:`32`)
  Maximum number of refinement iterations

:nml_n:`resolve_ctr` (default :nml_v:`.TRUE.`)
  Flag to resolve central evanescent region

:nml_n:`file` (default :nml_v:`''`)
   Name of file containing skeleton grid data (when
   :nml_n:`grid_source`\ =\ :nml_v:`'DETAIL_FILE'` or
   :nml_v:`'TEXT_FILE'`)

:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match

See the :ref:`spatial-grids` section for further details, in
particular a discussion of how the weighting (:nml_n:`w_*`) parameters
work.
