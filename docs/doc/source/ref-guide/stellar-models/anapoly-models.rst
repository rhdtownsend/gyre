.. _anapoly-models:

Analytic Polytropic Models
==========================

Setting the :nml_n:`model_type` parameter of the :nml_g:`model`
namelist group to :nml_v:`'ANAPOLY'` tells the frontend to create an
analytic polytropic stellar model. Because the structure of these model
can be computed analytically, there is no need to read from an
external file.

The :nml_n:`n_poly` parameter of the :nml_g:`model` namelist group
controls the polytropic index; valid choices are 0, 1 or 5. Likewise,
the :nml_n:`Gamma_1` parameter controls the first adiabatic index of
the model, while the :nml_n:`n`, :nml_n:`s` and :nml_n:`grid_type`
parameters control the model grid. See the :ref:`model-params` section
for further details.

It's possible to create truncated polytropic models, with a finite
surface pressure and density, by setting the :nml_n:`theta_s` to a
value other than zero (this parameter represents the desired surface
value of the polytropic dependent variable :math:`\theta`). A non-zero
:nml_n:`theta_s` is required when :nml_n:`n_poly`\ =\ :nml_v:`5`.
