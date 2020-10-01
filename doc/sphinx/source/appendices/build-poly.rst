.. _build-poly:

***********************************
Building Polytropic Structure Files
***********************************

This appendix describes the :program:`build_poly` executable, which
builds a composite polytropic stellar model and writes it to a file in
the :nml_v:`POLY` format (see the :ref:`poly-format` section).

Composite Polytropes
--------------------

A composite polytropic model comprises :math:`N` regions extending
from the origin out to the stellar surface. In the :math:`i`'th region
(:math:`1 \leq i \leq N`), the density and pressure are related by the
polytropic equation of state

.. math::

   P = K_{i} \rho^{(n_{i} + 1)/n_{i}}

where :math:`K_{i}` and the polytropic index :math:`n_{i}` are
constant across the region (but may differ from one region
to the next). At the :math:`N-1` boundaries between adjacent
regions, the pressure is required to be continuous but the
density may jump.

The structure of a composite polytrope is found by integrating the
composite Lane-Emden equation

.. math::

   \frac{1}{\xi^{2}} \deriv{}{z} \left( \xi^{2} \deriv{\theta_{i}}{\xi} \right) = - B_{i} \theta^{n_{i}}

in each region. Here, the independent variable is

.. math::

   \xi = \left[ \frac{K_{1} \rho_{\rm c}^{(1-n_{1})/n_{1}} (n_{1} + 1)}{4 \pi G} \right]^{-1/2} \, r,

where :math:`\rho_{\rm c}` is the central density of the star. In the
first region, the initial conditions are

.. math::

   \theta_{1}(0) = 1, \qquad \theta_{1}'(0) = 0

and :math:`B_{1} = 1`. In the subsequent regions, the initial conditions are

.. math::

   \theta_{i}(\xi_{i-1/2}) = 1, \qquad
   \theta'_{i}(\xi_{i-1/2}) = \frac{n_{i-1}+1}{n_{i}}
   \frac{\left[ \theta_{i}(\xi_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(\xi_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}}{t_{i-1}} \, \theta'_{i-1}(\xi_{i-1/2})

where :math:`\xi_{i-1/2}` is the coordinate of the boundary between the
:math:`i-1` and :math:`i` regions; moreover,

.. math::

   B_{i} = B_{i-1} \frac{n_{i-1}+1}{n_{i}+1}
   \frac{\left[ \theta_{i}(\xi_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(\xi_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}^{2}}{t_{i-1}^{2}}

and the recurrence for :math:`t_{i}` is

.. math::

   \ln t_{k} = \ln t_{i-1} + n_{i-1} \ln \theta_{i-1}(\xi_{i-1/2}) - n_{i} \ln \theta_{i}(\xi_{i-1/2}) + \Delta_{i-1/2},

with :math:`t_{1} = 1` and

.. math::

   \Delta_{i-1/2} = \ln \left[ \frac{\rho_{i}(\xi_{i-1/2})}{\rho_{i-1}(\xi_{i-1/2})}

quantifying the density jump at :math:`\xi_{i-1/2}`.

The surface of the composite polytropic model, :math:`\xi=\xi_{\rm
s}`, is defined implicitly by the boundary condition

.. math::

   \theta_{N}(\xi_{\rm s}) = 0.


Pre-Requisites
==============

In addition to GYRE's general pre-requisites (see the
:ref:`installation` chapter), :program:`build_poly` needs a
thread-safe version of the :netlib:`ODEPACK <odepack>` ordinary
differential integrator library. This library is shipped with version
20.3.2 (and more recent) of the `MESA SDK <mesa-sdk_>`__.

Compiling
=========

The :program:`build_poly` executable is automatically compiled when
GYRE is built, and installed in the :file:`{$GYRE_DIR}/bin` directory
(see the :ref:`installation` chapter).

Running
=======


