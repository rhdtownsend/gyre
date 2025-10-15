*************
Preliminaries
*************

Intended Audience
=================

This manual is aimed at a broad audience --- whether you're a GYRE
novice or a seasoned veteran, it provides you with the information
you'll need to get the most out of GYRE. However, it does presume some
experience with Unix command-line environments, and likewise some
basic familiarity with the subject of stellar oscillations. If you
need the former, then the Internet is your oyster; and for the latter,
we recommend the following online resources:

* Jørgen Christensen-Dalsgaard's `Lecture Notes on Stellar Oscillations <https://users-phys.au.dk/jcd/oscilnotes/Lecture_Notes_on_Stellar_Oscillations.pdf>`__;
* Gerald Handler's `Asteroseismology <https://arxiv.org/pdf/1205.6407.pdf>`__ article.

Obtaining GYRE
==============

The source code for GYRE is hosted in the :git:`rhdtownsend/gyre` git
repository on :git:`GitHub <>`. GYRE is free software: you can
redistribute it and/or modify it under the terms of the `GNU General
Public License <http://www.gnu.org/licenses/gpl-3.0.html>`__ as published
by the `Free Software Foundation <https://www.fsf.org/>`__, version 3.

.. _citing-gyre:

Citing GYRE
===========

If you use GYRE in your research, please cite one or more of the
relevant 'instrument' papers:

* :ads_citet:`townsend:2013` describes the basic operation of the code;
* :ads_citet:`townsend:2018` outlines capabilities for :ref:`non-adiabatic oscillation <non-ad-osc>` calculations;
* :ads_citet:`goldstein:2020` describes the :ref:`contour method <non-ad-contour>` for finding
  non-adiabatic modes;
* :ads_citet:`sun:2023` introduces the :ref:`gyre_tides frontend
  <frontends-gyre_tides>` for evaluating tidal responses.

If you find yourself using GYRE on a regular basis, you might consider
contributing to the project to ensure its long-term success. Options include

* contributing code to the project (e.g., via GitHub pull requests), to
  extend GYRE's capabilities;
* contributing documentation and tutorials to the project, to make GYRE more user-friendly;
* inviting the GYRE team to be co-authors on relevant papers;
* inviting the GYRE team to be co-investigators on relevant grant applications.

Development Team
================

GYRE remains under active development by the following team:

* `Rich Townsend <http://www.astro.wisc.edu/~townsend>`__ (University of Wisconsin-Madison); project leader
* `Warrick Ball <https://www.birmingham.ac.uk/staff/profiles/physics/ball-warrick.aspx>`__ (University of Birmingham)
* `Earl Bellinger <https://earlbellinger.com/>`__ (MPIA Garching)
* Zhao Guo (Cambridge University)
* Rianna Kuenzi (University of Wisconsin-Madison)
* Mathias Michielsen (KU-Leuven)
* Joel Ong (University of Hawaii-Manoa)
* `Jarret Rosenberg <https://www.physics.wisc.edu/directory/rosenberg-jarrett/>`__ (University of Wisconsin-Madison)
* Meng Sun (Northwestern University)
* Vincent Vanlaer (KU-Leuven)

Former developers include:

* Jacqueline Goldstein (MIT)
* Aaron Lopez

Also, the following people have made valuable contributions toward testing GYRE:

* Siemen Burssens (KU Leuven)
* Timothy Van Reeth (KU Leven)

Related Links
=============

* The `GYRE discussion forums <gyre-forums_>`__, the place to post
  feature requests and bug reports (don't send emails!).
* The `MESA Software Development Kit (SDK) <mesa-sdk_>`__, which
  provides the compilers and supporting libraries needed to build
  GYRE.
* The `MESA Stellar Evolution Code <mesa_>`__, which can generate
  stellar models readable by GYRE.

Acknowledgments
================

GYRE has been developed with financial support from the following grants:

* NSF awards AST-0908688, AST-0904607, ACI-1339606, ACI-1663696, and AST-1716436;
* NASA awards NNX14AB55G, NNX16AB97G, and 80NSSC20K0515.

GYRE has also benefited greatly from contributions (code, bug
reports, feature requests) from the academic community. Thanks, folks!
