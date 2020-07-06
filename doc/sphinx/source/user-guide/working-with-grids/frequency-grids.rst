.. _frequency-grids:

Frequency Grids
===============

GYRE searches for sign changes of the discriminant function
:math:`\Dfunc(\omega)` on a grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}` in the
dimensionless frequency :math:`\omega \equiv \sqrt{R^{3}/GM} \sigma`.
The computational cost of a calculation scales with the total
number of points :math:`M` in this grid, while the grid's resolution
--- i.e., the spacing between adjacent points --- impacts the
completeness of the modes found by GYRE (see the
:ref:`numerical-limits` section for a discussion of these behaviors in
the context of the stretched string BVP).

GYRE constructs a fresh frequency grid for each combination of
harmonic degree :math:`\ell` and azimuthal order :math:`m` specified
in the :nml_g:`mode` namelist groups (see the
:ref:`namelist-input-files` chapter for more details). 


     The starting
point for each of these grids is the *scaffold grid*, which comprises
the following:

* an inner point :math:`x=\xin`;
* an outer point :math:`x=\xout`;
* the subset of points of the input model grid satisfying :math:`\xin <
  x < \xout`

By default, :math:`\xin` and :math:`\xout` are obtained from the input
model grid as well, meaning that the scaffold grid is identical to the
model grid. However, either or both can be overridden using the
:nml_n:`x_i` and :nml_n:`x_o` parameters, respectively, of the
:nml_g:`grid` namelist group.



The grid is established 
