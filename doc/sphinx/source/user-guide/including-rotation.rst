.. _including-rotation:

******************
Including Rotation
******************

GYRE can partially include the Doppler and Coriolis effects arising
from stellar rotation, as described in detail in the
:ref:`rot-effects` section. Here, we briefly review the parameters
that influence how these effects are calculated.

Setting the Rotation Rate
-------------------------

There are a number of different ways to define the rotation angular
velocity :math:`\Omega`, via parameters in the :nml_g:`rot` namelist
group.

* If :nml_n:`Omega_rot_source` = :nml_v:`MODEL`, then :math:`\Omega`
  is obtained from the stellar model (for those :ref:`model formats
  <evol-models>` that include these data)

* If :nml_n:`Omega_rot_source` = :nml_v:`UNIFORM`, then uniform
  rotation is assumed, with the spatially constant :math:`\Omega` set
  by the :nml_n:`Omega_rot` and :nml_n:`Omega_rot_units` parameters

Enabling Doppler Effects
------------------------

Doppler effects are enabled by default whenever calculations are
performed with non-zero :math:`\Omega` and mode azimuthal order
:math:`m`. To disable these effects, set :nml_n:`m`\ = :nml_v:`0` in
the :nml_g:`mode` namelist group.

Enabling Coriolis Effects
-------------------------

Coriolis effects are disabled by default, but can be enabled using the
traditional approximation of rotation (TAR) by setting
:nml_n:`coriolis_method`\ =\ :nml_v:`'TAR'` in the :nml_g:`rot`
namelist group. The solution family is controlled by the
:nml_n:`rossby` parameter of the :nml_g:`rot` namelist group; set to
:nml_v:`.TRUE.` for the Rossby family, and to :nml_v:`.FALSE.` (the
default) for the gravito-acoustic family.



