.. _anapoly-models:

Analytic Polytropic Models
==========================

.. nml:group:: model
   :no-target:

Setting the :nml:option:`model_type` option of the :nml:group:`model`
namelist group to one of :nml:value:`'ANAPOLY_0'`,
:nml:value:`'ANAPOLY_1'`, or :nml:value:`'ANAPOLY_5'` tells the
frontend to create an analytic polytropic stellar model with the
indicated polytropic index (e.g., :nml:value:`ANAPOLY_1` has an index
:math:`\npoly=1`). Because the structure of these model can be
computed analytically, there is no need to read from an external
file. The :nml:value:`'ANAPOLY_5_1'` option is a special case,
constructed by matching an inner :math:`\npoly=5` region to an outer
:math:`\npoly=1` region; the location of the matching point is set by
the :nml:option:`x_match` option (see :ads_citealp:`eggleton:1998` for
details).

The :nml:option:`Gamma_1` option controls the first adiabatic index of
the model, while the :nml:option:`n`, :nml:option:`s` and :nml:option:`grid_type`
options control the model grid. See the :ref:`model-group` section
for further information.

It's possible to create truncated polytropic models, with a finite
surface pressure and density, by setting :nml:option:`theta_s` to a
value other than zero (this option represents the desired surface
value of the polytropic dependent variable :math:`\theta`). A non-zero
:nml:option:`theta_s` is required when using :nml:value:`'ANAPOLY_5'`
models, because the star would otherwise have an infinite radius.
