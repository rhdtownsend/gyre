.. _grid-group:

.. nml:group:: grid

Grid Namelist Group
===================

The :nml:group:`grid` namelist group controls how the spatial grid is
constructed. The following options are available:

.. nml:option:: scaffold_src
   :type: string
   :default: 'MODEL'

   Source for scaffold grid; one of:

   - :nml:value:`'MODEL'` : Obtained from the stellar model
   - :nml:value:`'FILE'` : Read from a file (see the :nml:option:`file` and :nml:option:`file_format` options)

.. nml:option:: x_i
   :type: real

   Inner boundary coordinate of calculation grid. If not specified,
   taken from the model grid

.. nml:option:: x_o
   :type: real

   Outer boundary coordinate of calculation grid. If not specified,
   taken from the model grid

.. nml:option:: w_osc
   :type: real
   :default: 0

   Oscillatory weight :math:`w_{\rm osc}`

.. nml:option:: w_exp
   :type: real
   :default: 0

   Exponential weight :math:`w_{\rm exp}`

.. nml:option:: w_ctr
   :type: real
   :default: 0

   Center weight :math:`w_{\rm ctr}`

.. nml:option:: w_thm
   :type: real
   :default: 0

   Thermal weight :math:`w_{\rm thm}`

.. nml:option:: w_str
   :type: real
   :default: 0

   Structural weight :math:`w_{\rm str}`

.. nml:option:: dx_min
   :type: real
   :default: SQRT(EPSILON(1._WP))

   Minimum spacing of grid points

.. nml:option:: dx_max
   :type: real
   :default: HUGE(0._WP)

   Maximum spacing of grid points

.. nml:option:: n_iter_max
   :type: integer
   :default: 32

   Maximum number of refinement iterations

.. nml:option:: resolve_ctr
   :type: logical
   :default: .TRUE.

   Flag to resolve central evanescent region

.. nml:option:: file
   :type: string
   :default: ''

   Name of file containing scaffold grid data. Used only when
   :nml:option:`scaffold_src` = :nml:value:`'FILE'`

.. nml:option:: file_format
   :type: string
   :default: ''

   Format of file containing scaffold grid data; one of

   - :nml:value:`'TEXT'`: text file with one abscissa value per line
   - :nml:value:`'DETAIL'`: :ref:`detail <detail-files>` file with
     abscissa values provided in ``x`` dataset

   Used only when :nml:option:`scaffold_src` = :nml:value:`'FILE'`

.. nml:option:: tag_list
   :type: string
   :default: ''

   Comma-separated list of :nml:option:`tag <mode.tag>` values to
   match; matches all if left blank

See the :ref:`spatial-grids` section for further details, in
particular a discussion of how the weight options
(:nml:option:`w_osc`, :nml:option:`w_exp`, etc.) work.
