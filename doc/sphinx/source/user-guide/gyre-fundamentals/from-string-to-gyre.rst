.. _from-string-to-gyre:

From Stretched String to GYRE
=============================

The stretched-string BVP explored in this chapter provides a very nice
analog of the numerical technique that GYRE uses to solve the
oscillation equations. In this section, let's review where this

Similarities in GYRE
--------------------

* Separation of variables
* Discretization
* Root bracketing for adiabatic
* Root refinement
* Failure scenarios

Differences in GYRE
-------------------

* More variables
* Non-uniform spatial and frequency grids
* The stretched string is discretized using a second-order form of the wave equation. GYRE breaks this up into 1st-order forms
* The discretizaton scheme is more complex
* Complex root solving for non-adiabatic
  
