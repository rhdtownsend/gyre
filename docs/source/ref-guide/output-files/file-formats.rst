.. _file-formats:

File Formats
============

The format of summary and detail files depends on the value of the
:nml_n:`summary_file_format` and :nml_n:`detail_file_format` parameters
in the :nml_g:`ad_output` and :nml_g:`nad_output` namelist groups (see
the :ref:`output-params` section). Possible choices are:

* :nml_v:`'HDF'` : A binary format based on the `HDF5
  <https://support.hdfgroup.org/HDF5/whatishdf5.html>`__ format

* :nml_v:`'TXT'` : A text format modeled after
  MESA's `profile file format <http://mesa.sourceforge.net/output.html>`__

For both formats, the data stored in the files come in two flavors ---
scalars (a single value) and arrays (a sequence of values). Files in
either format can be read in Python using the
:external:py:func:`pygyre.read_output` function from :git:`PyGYRE
<rhdtownsend/pygyre>` (see the :ref:`interpreting-output` chapter for
examples).

.. _hdf-format:

HDF Format
----------

HDF-format output files adhere to the following conventions:

* All data objects are attached to the root HDF5 group (`/`)
* Attributes are used to store scalar data
* Datasets are used to store array data
* Real values are written with type `H5T_IEEE_F64LE` when GYRE is
  compiled in double precision (the default), and type
  `H5T_IEEE_F32LE` otherwise
* Integer values are written with type `H5T_STD_I32LE`
* Complex values are written as a compound type, composed of a real
  component `re` and an imaginary component `im`; the types of
  these components are the same as for real values

.. _txt-format:

TXT Format
----------

TXT-format files adhere to the following conventions:

* The first three lines contain the scalar data:

  * The first line contains the column numbers for the scalar data,
    starting at 1
  * The second line contains the column labels for the scalar data
  * The third line contains the actual scalar data values

* The subsequent lines contain the array data:

  * The fourth line contains the column numbers for the array data,
    starting at 1
  * The fifth line contains the the column labels for the array data
  * The sixth and subsequent lines contain the actual array data (one
    line per array element)

* Complex values are written as two columns, with the first column
  containing the real component and the second the imaginary component
