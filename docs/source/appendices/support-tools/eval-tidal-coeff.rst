.. _support-tools-eval-tidal-coeff:

eval_tidal_coeff
================

.. program:: eval_tidal_coeff

The :program:`eval_tidal_coeff` tool evaluates one of a number of
coefficients used in tidal expansions, as follows:

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Label
     - Symbol
     - Description
   * - :code:`cbar`
     - :math:`\cbar_{\ell,m,k}`
     - expansion coefficient for tidal potential
   * - :code:`Gbar_1`
     - :math:`\Gbar^{(1)}_{\ell,m,k}`
     - expansion coefficient for secular rate-of-change of longitude of periapsis
   * - :code:`Gbar_2`
     - :math:`\Gbar^{(2)}_{\ell,m,k}`
     - expansion coefficient for secular rate-of-change of semi-major axis
   * - :code:`Gbar_3`
     - :math:`\Gbar^{(3)}_{\ell,m,k}`
     - expansion coefficient for secular rate-of-change of eccentricity
   * - :code:`Gbar_4`
     - :math:`\Gbar^{(4)}_{\ell,m,k}`
     - expansion coefficient for secular rate-of-change of spin
   * - :code:`X`
     - :math:`X^{-(\ell+1),-m}_{-k}`
     - Hansen coefficient
   * - :code:`Y`
     - :math:`\real[Y^{m}_{\ell}(\pi/2,0)]`
     - spherical harmonic at equator/prime meridian
   * - :code:`Y*`
     - :math:`\real[Y^{m*}_{\ell}(\pi/2,0)]`
     - complex conjugate spherical harmonic at equator/prime meridian

The full definitions of these coefficients can be found in :ads_citet:`sun:2023`.

Syntax
------

:program:`eval_tidal_coeff` :option:`coeff` :option:`R_a` :option:`e` :option:`l` :option:`m` option:`k`

Parameters
----------

.. option:: coeff

   Coefficient to calculate. For possible choices, see the label field in the table above.

.. option:: R_a

   Ratio of stellar radius :math:`R` to orbital semi-major axis :math:`a`

.. option:: e

   Orbital eccentricity :math:`e`

.. option:: l

   Harmonic degree :math:`\ell`

.. option:: m

   Azimuthal order :math:`m`

.. option:: k

   Fourier index :math:`k`

Output
------

The coefficient is printed to standard output.

