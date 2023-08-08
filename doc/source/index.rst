#################################
The GYRE Stellar Oscillation Code
#################################

GYRE is a *stellar oscillation code*. Given an input stellar model,
GYRE calculates the eigenfrequencies and eigenfunctions for the normal
oscillation modes of the model. These data can be put to a variety of
uses; the most common is to compare them against observed oscillation
frequencies of a star, allowing constraints on the star's fundamental
parameters (mass, radius, etc.)  to be established --- the discipline
of *asteroseismology*.

GYRE also supports other, related kinds of calculation. One example is
evaluating the response of a star to tidal disturbances produced by an
orbiting companion; because this is an instance of *forced* stellar
oscillations, similar numerical techniques can be brought to bear.

About this Manual
=================

This manual is divided into two main parts: the :ref:`user-guide`,
which covers basic installation and usage of GYRE, and the
:ref:`ref-guide`, which documents all of GYRE's options in
detail. Supplementary material can be found in the :ref:`appendices`.

.. toctree::
   :caption: User Guide
   :name: user-guide
   :maxdepth: 2

   user-guide/preliminaries.rst
   user-guide/quick-start.rst
   user-guide/example-walkthrough.rst
   user-guide/frontends.rst
   user-guide/numerical-methods.rst
   user-guide/interpreting-output.rst
   user-guide/understanding-grids.rst
   user-guide/working-with-tags.rst
   user-guide/advanced-usage.rst
   user-guide/faq.rst
   user-guide/troubleshooting.rst
   
.. toctree::
   :caption: Reference Guide
   :name: ref-guide
   :maxdepth: 2

   ref-guide/installation.rst
   ref-guide/input-files.rst
   ref-guide/output-files.rst
   ref-guide/stellar-models.rst
   ref-guide/osc-equations.rst

.. toctree::
   :caption: Appendices
   :name: appendices

   appendices/comp-ptrope.rst
   appendices/build-poly.rst
   appendices/eval-lambda.rst
   
