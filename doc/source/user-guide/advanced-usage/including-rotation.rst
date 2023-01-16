.. _rotation:

Including Rotation
==================

This section discusses how to include the effects of rotation in GYRE
calculations. Further details of GYRE's treatment of rotation are
presented in the :ref:`osc-eqns` and :ref:`tidal-eqns` chapters.

Setting the Rotation Rate
-------------------------

There are two different ways to define the rotation angular velocity
:math:`\Omega`, via parameters in the :nml_g:`rot` namelist group.

* If :nml_n:`Omega_rot_source` = :nml_v:`MODEL`, then differential
  rotation is assumed with a spatially varying :math:`\Omega`
  obtained from the stellar model. If the model doesn't have this
  capability (see the :ref:`model-caps` section), then :math:`\Omega`
  is set to zero throughout the star.
  
* If :nml_n:`Omega_rot_source` = :nml_v:`UNIFORM`, then uniform
  rotation is assumed with a spatially constant :math:`\Omega` set
  by the :nml_n:`Omega_rot` and :nml_n:`Omega_rot_units` parameters.

Enabling Doppler Effects
------------------------

:ref:`Doppler effects <osc-rot-doppler>` are enabled by default
whenever calculations are performed with non-zero :math:`\Omega` and
mode azimuthal order :math:`m`. To disable these effects, set
:nml_n:`m`\ = :nml_v:`0` in the :nml_g:`mode` namelist group.

Enabling Coriolis Effects
-------------------------

:ref:`Coriolis effects <osc-rot-coriolis>` are disabled by default,
but can be enabled using the traditional approximation of rotation
(TAR) by setting :nml_n:`coriolis_method`\ =\ :nml_v:`'TAR'` in the
:nml_g:`rot` namelist group. The solution family is controlled by the
:nml_n:`rossby` parameter of the :nml_g:`rot` namelist group; set to
:nml_v:`.TRUE.` for the Rossby family, and to :nml_v:`.FALSE.` (the
default) for the gravito-acoustic family.
