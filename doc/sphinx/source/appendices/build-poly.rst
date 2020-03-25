.. _build-poly:

***********************************
Building Polytropic Structure Files
***********************************

This appendix describes the :program:`build_poly` executable, which creates polytropic structure
files for use with GYRE.

Pre-Requisites
==============

In addition to GYRE's general pre-requisites (see the
:ref:`installation` chapter), :program:`build_poly` needs a
thread-safe version of the :netlib:`ODEPACK <odepack>` ordinary
differential integrator library. This library is shipped with version
20.3.2 (and more recent) of the MESA SDK.


On Linux and MacOS platforms, this library is bundled together with all of the other pre-requisites in the Madison Software Development Kit (SDK), which can be downloaded from
the `MESA SDK <mesa-sdk_>`__ homepage. Using this SDK is strongly
recommended.


It does this by solving the Lane-Emden
equation

.. math::

   \frac{1}{\xi^{2}} \deriv{}{\xi} \left( \xi^{2} \deriv{\Theta}{\xi} \right) = - \Theta^{\npoly}

for constant polytropic index :math:`\npoly`. Here, the dependent
variable :math:`\Theta` represents the local density :math:`\rho` via

.. math::

   \Theta = \left( \frac{\rho}{\rho_{\rm c} \right)^{1/\npoly},

where :math:`\rho_{\rm c}` is the central density; and the independent
variable :math:\`Theta` is related to the radial coordinate :math:`r`
via

.. math::

   \frac{r}{R} = \frac{\xi}{\xi_{1}},

where :math:`R` and :math:`\xi_{1}` are the surface values of
:math:`r` and :math:`\xi`, respectively. The Lane-Emden is solved
subject to the boundary conditions

.. math::

   \theta(0) = 1, \qquad \theta(\xi_{1}) = 0

Compile
-------

