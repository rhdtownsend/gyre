.. _support-tools-eval-tidal-coeff:

eval_tidal_coeff
================

.. program:: eval_tidal_coeff

Synopsis
--------

.. code-block:: text

   eval_lambda ( cbar | Gbar_{1..4} | X | Y | Y^* ) [options]

Description
-----------

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
     - :math:`X^{n,m}_{k}`
     - Hansen coefficient
   * - :code:`Y`
     - :math:`\real[Y^{m}_{\ell}(\pi/2,0)]`
     - spherical harmonic at equator/prime meridian
   * - :code:`Y^*`
     - :math:`\real[Y^{m*}_{\ell}(\pi/2,0)]`
     - complex conjugate spherical harmonic at equator/prime meridian

The coefficient is written to standard output.

Options
-------

.. option:: -h, --help

   Print a summary of options.

.. option:: --R/a=R_A

   Ratio of stellar radius :math:`R` to orbital semi-major axis :math:`a`.

.. option:: -e, --e=E

   Orbital eccentricity :math:`e`.

.. option:: -n, --n=N

   Hansen degree :math:`n`.

.. option:: -l, --l=L

   Harmonic degree :math:`\ell`.

.. option:: -m, --m=M

   Azimuthal order :math:`m`

.. option:: -k, --k=K

   Fourier index :math:`k`.

Output
------

The coefficient is printed to standard output.

