.. _from-string-to-gyre:

From Stretched String to GYRE
=============================

The numerical technique used in this chapter to solve the
stretched-string BVP provides a strong analog to how GYRE solves the
oscillation equations. The full, nasty details of GYRE's approach are
laid out in XXXX

In this section, let's review the similarities
and differences between the two.

Similarities
------------

GYRE shares the following similarities to 

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
  
