.. _anapoly-models:

Analytic Polytropic Models
==========================

Setting the :nml_n:`model_type` parameter of the :nml_g:`model`
namelist group to one of :nml_v:`'ANAPOLY_0'`, :nml_v:`'ANAPOLY_1'`,
or :nml_v:`'ANAPOLY_5'` tells the frontend to create an analytic
polytropic stellar model with the indicated polytropic index (e.g.,
:nml_v:`ANAPOLY_1` has an index :math:`\npoly=1`). Because the
structure of these model can be computed analytically, there is no
need to read from an external file. The :nml_v:`'ANAPOLY_5_1'` option
is a special case, constructed by matching an inner :math:`\npoly=5`
region to an outer :math:`\npoly=1` region.

The :nml_n:`Gamma_1` parameter controls the first adiabatic index of
the model, while the :nml_n:`n`, :nml_n:`s` and :nml_n:`grid_type`
parameters control the model grid. See the :ref:`model-params` section
for further details.

It's possible to create truncated polytropic models, with a finite
surface pressure and density, by setting the :nml_n:`theta_s` to a
value other than zero (this parameter represents the desired surface
value of the polytropic dependent variable :math:`\theta`). A non-zero
:nml_n:`theta_s` is required when using :nml_v:`'ANAPOLY_5'` models,
because the star would otherwise have an infinite radius.
