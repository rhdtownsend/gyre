.. _tidal-rot:

Rotation Effects
================

Including rotation in tidal response calculations is complicated by
the fact that the :ref:`traditional approximation of rotation
<osc-rot-tar>` doesn't play well with forcing by the tidal potential
:math:`\PhiT`. As a result, for now only the :ref:`Doppler shift
<osc-rot-doppler>` due to rotation is included in
:program:`gyre_tides` calculations. This means that the frequency
:math:`\sigma` appearing in the :ref:`tidal equations
<tidal-sep-eqns>` is replaced by the co-rotating forcing frequency

.. math::
   :label: e:tidal-sigmac

   \sigmac \equiv k \Oorb - m \Omega.
