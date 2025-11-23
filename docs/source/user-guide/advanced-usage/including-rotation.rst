.. _rotation:

Including Rotation
==================

This section discusses how to include the effects of rotation in GYRE
calculations. See the :ref:`osc-rot` section for further details of
GYRE's rotation treatment.

Setting the Rotation Rate
-------------------------

There are two different ways to define the rotation angular frequency
:math:`\Orot`, via parameters in the :nml_g:`rot` namelist group.

* If :nml_n:`Omega_rot_source` = :nml_v:`'MODEL'`, then differential
  rotation is assumed with a spatially varying :math:`\Orot`
  obtained from the stellar model. If the model doesn't have this
  capability (see the :ref:`model-caps` section), then :math:`\Orot`
  is set to zero throughout the star.

* If :nml_n:`Omega_rot_source` = :nml_v:`'UNIFORM'`, then uniform
  rotation is assumed with a spatially constant :math:`\Orot` set by
  the :nml_n:`Omega_rot` and :nml_n:`Omega_rot_units` parameters.

Incorporating Doppler Effects
-----------------------------

The :ref:`osc-rot-doppler` effect is incorporated automatically
whenever calculations are performed with non-zero :math:`\Orot` and
mode azimuthal order :math:`m`. To disable this effect, set
:nml_n:`m`\ = :nml_v:`0` in the :nml_g:`mode` namelist group.

Incorporating Coriolis Effects
------------------------------

Incorporating the effects of the Coriolis force can be done using a
:ref:`perturbative treatment <osc-rot-coriolis-p>` or a
:ref:`non-perturbative treatment <osc-rot-coriolis-np>`. In the former
case the effects are be applied as a post-calculation correction to
non-rotating eigenfrequencies (see the :nml_v:`domega_rot` output item
in the :ref:`summary-files` and :ref:`detail-files` sections). In the
latter case, the traditional approximation of rotation (TAR) can be
enabled by setting the :nml_n:`lambda_method` in the :nml_g:`osc`
namelist group, to :nml_v:`'TAR-GRAVITY'` (for gravito-acoustic modes)
or to :nml_v:`'TAR-ROSSBY'` (for Rossby modes).
