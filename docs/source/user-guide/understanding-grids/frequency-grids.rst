.. _freq-grids:

Frequency Grids
===============

The :program:`gyre` frontend evaluates its discriminant function
:math:`\Dfunc(\omega)` on a grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}` in the
dimensionless frequency, and scans for changes in the sign of
:math:`\Dfunc(\omega)` that are indicative of a bracketed root.  The
computational cost of a calculation scales with the total number of
points :math:`M` in this grid, while the grid's resolution --- i.e.,
the spacing between adjacent points --- impacts the completeness of
the modes found by :program:`gyre`. (See the :ref:`numerical-limits`
section for a discussion of these behaviors).

A fresh frequency grid is constructed for each iteration of the main
computation loop (see the flow-chart in the :ref:`gyre
<frontends-gyre>` section). This is done under the control of the
:nml_g:`scan` namelist groups; there must be at least one of these,
subject to the tag matching rules (see the :ref:`working-with-tags`
chapter). Each :nml_g:`scan` group creates a separate frequency grid;
these are then combined and the discriminant function is evaluated on
the merged grid.

Grid Types
----------

The :nml_n:`grid_type` parameter of the :nml_g:`scan` namelist group
controls the overall distribution of points in a frequency grid. There
are currently three options:

.. _linear-freq-grid:

Linear Grid
~~~~~~~~~~~

When :nml_nv:`grid_type = 'LINEAR'`, :program:`gyre` first evaluates a
sequence of dimensionless angular frequencies in the grid reference
frame according to the formula

.. math::

   \omega^{\rm g}_{i} = \frac{1}{M-1} \left[ (M - i)\, \omega^{\rm g}_{\rm min}  + (i - 1) \, \omega^{\rm g}_{\rm max} \right]
   \qquad i = 1,2,\ldots,M.

(here, the superscript 'g' indicates that these are frequencies in the
grid reference frame). Then, it transforms from the grid frame to the
inertial reference frame via

.. math::

   \omega_{i} = \omega^{\rm g}_{i} + \Delta \omega

where :math:`\Delta\omega` is the frequency shift that transforms from
the grid frame to the inertial frame. The actual value of this shift
depends on the :nml_n:`grid_frame` parameter:

* When :nml_n:`grid_frame = 'INERTIAL'`, the shift is :math:`\Delta
  \omega = 0`; the grid frame and the inertial frame coincide.

* When :nml_n:`grid_frame = 'COROT_I'`, the shift is :math:`\Delta
  \omega = m \, \Orot^{\rm i} \sqrt{R^{3}/GM}`, where
  :math:`\Orot^{\rm i}` is the rotation angular frequency at the
  inner boundary of the spatial grid; the grid frame coincides with
  the local co-rotating frame at that boundary.

* When :nml_n:`grid_frame = 'COROT_O'`, the shift is :math:`\Delta
  \omega = m \, \Orot^{\rm o} \sqrt{R^{3}/GM}`, where
  :math:`\Orot^{\rm o}` is the rotation angular frequency at the outer
  boundary of the spatial grid; the grid frame coincides with the
  local co-rotating frame at that boundary.

The range spanned by the frequency grid, in the grid frame, is set by
:math:`\omega^{\rm g}_{\rm min}` and :math:`\omega^{\rm g}_{\rm max}`. These are
evaluated via

.. math::

   \omega^{\rm g}_{\rm  min} = \frac{f_{\rm min}}{\widehat{f}_{\rm min}} + \delta \omega - \Delta \omega,
   \qquad \qquad
   \omega^{\rm g}_{\rm max} = \frac{f_{\rm max}}{\widehat{f}_{\rm max}} + \delta \omega - \Delta \omega,

where :math:`f_{\rm min,max}` are user-definable,
:math:`\widehat{f}_{\rm min,max}` will be discussed below in the
:ref:`freq-units` section, and :math:`\delta\omega` is the frequency
shift that transforms from the frame in which :math:`f_{\rm min,max}`
are defined to the inertial frame. The actual value of this shift depends
on the :nml_n:`freq_frame` parameter, which behaves analogously to the
:nml_n:`grid_frame` parameter discussed above.

.. _inverse-freq-grid:

Inverse Grid
~~~~~~~~~~~~

When :nml_nv:`grid_type = 'INVERSE'`, :program:`gyre` first evaluates a sequence
of dimensionless angular frequencies in the grid reference frame
according to the formula

.. math::

   \omega^{\rm g}_{i} = (M-1) \left[ \frac{(M - i)}{\omega^{\rm g}_{\rm min}}  + \frac{(i - 1)}{\omega^{\rm g}_{\rm max}} \right]^{-1}
   \qquad i = 1,2,\ldots,M.

The grid creation then proceeds as described above in the :ref:`linear-freq-grid` section.

File Grid
~~~~~~~~~

When :nml_nv:`grid_type = 'FILE'`, :program:`gyre` first reads a
sequence of dimensioned frequencies
:math:`\{f_{1},f_{2},\ldots,f_{M}\}` from an external file named by
the :nml_n:`grid_file` parameter. This file is a single-column ASCII
text table; the number of points :math:`M` is determined implicitly
from the number of lines in the file. Then, it transforms these
frequencies via

.. math::

   \omega_{i} = \frac{f_{j}}{\widehat{f}} + \delta \omega,

where :math:`\widehat{f}` will be discussed below in the
:ref:`freq-units` section, and :math:`\delta\omega` is the frequency
shift that transforms from the frame in which :math:`f` is defined to
the inertial frame. The actual value of this shift depends on the
:nml_n:`freq_frame` parameter, which behaves analogously to the
:nml_n:`grid_frame` parameter discussed above.

.. _freq-units:

Frequency Units
---------------

In the expressions above, terms of the form :math:`f/\widehat{f}` are used
to transform a dimensioned frequency :math:`f` into a dimensionless
one :math:`\omega`. The scale factor :math:`\widehat{f}` depends on the
:nml_n:`freq_units` parameter. Thus, for example, if
:nml_nv:`freq_units = 'UHZ'`, then :math:`f` is treated as a linear
frequency expressed in :math:`{\rm \mu Hz}`, and the scale factor is set by

.. math::

   \widehat{f} = \sqrt{\frac{GM}{R^{3}}} \frac{1}{2\pi\,{\rm \mu Hz}}

(the factor of :math:`2\pi` comes from the transformation between linear
and angular frequency).

The full set of values supported by the :nml_n:`freq_units` parameter
is listed in the :ref:`scan-params` section.

Namelist Parameters
-------------------

The full set of parameters supported by the :nml_g:`scan` namelist
group is listed in the :ref:`scan-params` section. However, the table
below summarizes the mapping between the user-definable controls
appearing in the expressions above, and the corresponding namelist
parameters:

.. list-table::
   :widths: 30 30
   :header-rows: 1

   * - Symbol
     - Parameter
   * - :math:`f_{\rm min}`
     - :nml_n:`freq_min`
   * - :math:`f_{\rm max}`
     - :nml_n:`freq_max`
   * - :math:`M`
     - :nml_n:`n_freq`

Recommended Values
------------------

The default values :nml_nv:`freq_min=1`, :nml_nv:`freq_max=10`,
:nml_nv:`n_freq=10`, together with :nml_nv:`grid_type='LINEAR'` are
sufficient to find *some* modes --- although unlikely the modes that
you want. Choosing good values for these parameters requires some
degree of judgment, but here are some suggestions:

* The number of points in the frequency grid should be a factor of
  2--3 larger than the number of modes you expect :program:`gyre` will
  find. This is to ensure that the frequency spacing of the grid is
  everywhere smaller than the anticipated eigenfrequency spacing
  between adjacent modes (see the :ref:`numerical-limits` section for
  further discussion).

* The distribution of points in the frequency grid should follow
  anticipated distribution of mode frequencies; this again is to
  ensure adequate frequency resolution. For p modes, which tend toward
  a uniform frequency spacing in the asymptotic limit of large radial
  order, you should chose :nml_nv:`grid_type = 'LINEAR'`;
  likewise, for g modes, which tend toward a uniform period spacing in
  the asymptotic limit, you should choose :nml_nv:`grid_type = 'INVERSE'`.

* When modeling rotating stars, you should choose :nml_nv:`grid_frame
  = 'COROT_I'` or :nml_nv:`grid_frame = 'COROT_O'`, because the
  asymptotic behaviors mentioned above apply in the co-rotating
  reference frame rather than the inertial one.


