.. _working-with-grids:

******************
Working with Grids
******************

Spatial Grids
=============

As discussed in the :ref:`gyre-fundamentals` chapeter, GYRE discretizes
the oscillation equations on a spatial grid. This grid spans a range
:math:`[x_{\rm i},x_{\rm o}]` in the dimensionless radial coordinate
:math:`x \equiv r/R`. The grid resolution 

GYRE solves the pulsation equations on a spatial grid --- a set of discrete points spanning some range of values [xi,xo] in the dimensionless radial coordinate x=r/R∗. This grid must be fine enough to adequately resolve the wavefunctions of the pulsation modes being sought; but this requirement must be balanced by the fact that the computation time scales (approximately linearly) with the number of points in the grid.

Frequency Grids
===============
