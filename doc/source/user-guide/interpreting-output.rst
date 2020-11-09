.. _interpreting-output:

*************************
Interpreting Output Files
*************************

This chapter reviews the summary and detail output files written
during a GYRE run, and demonstrates how to read and plot them in
`Python <https://www.python.org>`__. Further information about these
files is provided in the :ref:`output-files` chapter.

File Categories
===============

Summary files collect together global properties, such as
eigenfrequencies and radial orders, of all modes found. By contrast, a
detail file stores spatial quantities, such as eigenfunctions and
differential inertias, for an individual mode. The choice of which
specific data actually appear in output files isn't hardwired, but
rather determined by the :nml_n:`summary_item_list` and
:nml_n:`mode_item_list` parameters of the :nml_g:`ad_output` and
:nml_g:`nad_output` namelist groups. Changing these parameters allows
you to tailor the files to contain exactly the data you need.

File Formats
============

Summary and detail files are written by GYRE in either TXT or HDF
format. Files in the TXT format are human-readable, and can be
reviewed on-screen or in a text editor; whereas files in the HDF
format are intended to be accessed through a suitable `HDF5
<https://support.hdfgroup.org/HDF5/whatishdf5.html>`__
interface. Unless there's a good reason to use TXT format, HDF format
is preferable; it's portable between different platforms, and takes up
significantly less disk space

PyGYRE
======

:git:`PyGYRE <rhdtownsend/pygyre>` is a Python package, maintained
separately from GYRE, providing a set of routines that greatly
simplify the analysis of summary and detail files. Detailed
information about PyGYRE can be found in the `full documentation
<https://pygyre.readthedocs.io/en/latest/>`__; here, we demonstrate
how to use it to read and plot the output files from the
:ref:`walkthrough` section.

As a preliminary step, you'll need to install PyGYRE from the `Python
Package Index (PyPI) <https://pypi.org/>`__. This can be done using
the :command:`pip` command, via

.. substitution-prompt:: bash

   pip install pygyre

If PyGYRE is already installed, you can upgrade to a more-recent
version via

.. substitution-prompt:: bash

   pip install --upgrade pygyre

Analyzing a Summary File
========================

With PyGYRE installed, change into your :ref:`work directory
<walkthrough-work>` and fire up your preferred interactive Python
environment (e.g., `Jupyter <https://jupyter.org/>`__). Import PyGYRE
and the other modules needed for plotting:

.. code::

  # Import modules

  import pygyre as pg
  import matplotlib.pyplot as plt
  import numpy as np

(you may want to directly cut and paste this code). Next, read the
summary file in the work directory into the variable `s`:

.. code::

   # Read data from a GYRE summary file

   s = pg.read_output('summary.h5')

The ``pg.read_output`` function is able to read both TXT- and HDF-format
files, returning the data in a ``Table`` object (from the `Astropy
<https://www.astropy.org/>`__ project). To inspecting the data on-screen,
simply evaluate the table:

.. code::

   # Inspect the data

   s

From this, you'll see that there are three columns in the table,
containing the harmonic degree ``l``, radial order ``n_pg`` and
frequency ``freq`` of each mode found during the GYRE run.

Next, plot the frequencies against radial orders via

.. code::

   # Plot the data

   plt.figure()

   plt.plot(s['n_pg'], s['freq'].real)

   plt.xlabel('n_pg')
   plt.ylabel('Frequency (cyc/day)')

(the values in the ``freq`` column are complex, and we plot the real
part). The plot should look something along the lines of
:numref:`fig-freq`.

.. _fig-freq:

.. figure:: interpreting-output/fig_freq.svg
   :alt: Plot showing mode frequency versus radial order
   :align: center

   (:download:`Source <interpreting-output/fig_freq.py>`)

The straight line connecting the two curves occurs because we are
plotting both the dipole and quadrupole modes together. To separate
them, the table rows can be grouped by harmonic degree:

.. code::

   # Plot the data, grouped by harmonic degree

   plt.figure()

   sg = s.group_by('l')

   plt.plot(sg.groups[0]['n_pg'], sg.groups[0]['freq'].real, label=r'l=1')
   plt.plot(sg.groups[1]['n_pg'], sg.groups[1]['freq'].real, label=r'l=2')

   plt.xlabel('n_pg')
   plt.ylabel('Frequency (cyc/day)')

   plt.legend()

The resulting plot, in :numref:`fig-freq-grouped` looks much better.
   
.. _fig-freq-grouped:

.. figure:: interpreting-output/fig_freq_grouped.svg
   :alt: Plot showing mode frequency versus radial order
   :align: center

   (:download:`Source <interpreting-output/fig_freq_grouped.py>`)

Analyzing a Detail File
=======================

Now let's take a look at one of the detail files, for the mode with
:math:`\ell=1` and :math:`n_{\rm pg}=-7`. As with the summary file,
``pg.read_output`` can be used to read the file data into a ``Table``
object:

.. code::
   
   # Read data from a GYRE detail file

   d = pg.read_output('detail.l1.n-7.h5')

Inspecting the data using

.. code::

   # Inspect the data

   d

shows there are 7 columns: the fractional radius ``x``, the radial
displacement eigenfunction ``xi_r``, the horizontal displacement
eigenfunction ``xi_h``, and 4 further columns storing structure
coefficients (see the :ref:`detail-files` section for descriptions of
these data). Plot the two eigenfunctions using the code

.. code::

   # Plot displacement eigenfunctions

   plt.figure()

   plt.plot(d['x'], d['xi_r'].real, label='xi_r')
   plt.plot(d['x'], d['xi_h'].real, label='xi_h')

   plt.xlabel('x')

   plt.legend()

.. _fig-disp-eigfunc:

.. figure:: interpreting-output/fig_disp_eigfunc.svg
   :alt: Plot showing displacement eigenfunctions versus fractional radius
   :align: center

   The radial (:math:`\txir`) and horizontal (:math:`\txih`)
   displacement eigenfunctions of the :math:`\ell=1`, :math:`n_{\rm
   pg}=-7` mode, plotted against the fractional radius :math:`x`.
   (:download:`Source <interpreting-output/fig_disp_eigfunc.py>`)

The plot should look something along the lines of
:numref:`fig-disp-eigfunc`. From this figure , we see that the radial
wavelengths of the eigenfunctions become very short around a
fractional radius :math:`x \approx 0.125`. To figure out why this is, we
can take a look at the star's propagation diagram:

.. code::

   # Evaluate characteristic frequencies

   l = d.meta['l']
   omega = d.meta['omega']

   x = d['x']
   V = d['V_2']*d['x']**2
   As = d['As']
   c_1 = d['c_1']
   Gamma_1 = d['Gamma_1']

   d['N2'] = d['As']/d['c_1']
   d['Sl2'] = l*(l+1)*Gamma_1/(V*c_1)

   # Plot the propagation diagram

   plt.figure()

   plt.plot(d['x'], d['N2'], label='N^2')
   plt.plot(d['x'], d['Sl2'], label='S_l^2')

   plt.axhline(omega.real**2, dashes=(4,2))

   plt.xlabel('x')
   plt.ylabel('omega^2')

   plt.ylim(5e-2, 5e2)
   plt.yscale('log')

Note how we access the mode harmonic degree ``l`` and dimensionless
eigenfrequency ``omega`` through the table metadata dict
``d.meta``. The resulting plot (cf. :numref:`fig-prop-diag`) reveals
that the Brunt-Väisälä frequency squared is large around :math:`x
\approx 0.125`; this feature is a consequence of the molecular weight
gradient zone outside the star's convective core, and results in the
short radial wavelengths seen there in :numref:`fig-disp-eigfunc`.

.. _fig-prop-diag:

.. figure:: interpreting-output/fig_prop_diag.svg
   :alt: Plot showing propagation diagram
   :align: center

   Propagation diagram for the :math:`5\,\Msun` model, plotting the
   squares of the Brunt-Väisälä (:math:`N^{2}`) and Lamb
   (:math:`S_{\ell}^{2}`) frequencies versus fractional radius
   :math:`x`. The horizontal dashed line shows the frequency squared
   :math:`\omega^{2}` of the :math:`\ell=1`, :math:`n_{\rm pg}=-7`
   mode shown in :numref:`fig-disp-eigfunc`. Regions where
   :math:`\omega^{2}` is smaller (greater) than both :math:`N^{2}` and
   :math:`S_{\ell}^{2}` are gravity (acoustic) propagation regions;
   other regions are evanescent. (:download:`Source
   <interpreting-output/fig_prop_diag.py>`)
   


