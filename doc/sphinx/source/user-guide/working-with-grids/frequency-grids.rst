.. _frequency-grids:

Frequency Grids
===============

GYRE searches for sign changes of the discriminant function
:math:`\disc(\omega)` on a grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}` in the
dimensionless frequency :math:`\omega \equiv \sqrt{R^{3}/GM} \sigma`.
The computational cost of a calculation scales with the total
number of points :math:`M` in this grid, while the grid's resolution
--- i.e., the spacing between adjacent points --- impacts the
completeness of the modes found by GYRE (see the
:ref:`numerical-limits` section for a discussion of these behaviors in
the context of the stretched string BVP).

The grid 
