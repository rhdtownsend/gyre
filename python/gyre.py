# Module   : gyre
# Purpose  : process GYRE output
#
# Copyright 2013 Rich Townsend
#
# This file is part of GYRE. GYRE is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, version 3.
#
# GYRE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Imports

import h5py
import numpy as np

# Read a GYRE HDF5-format output file

def read_output (filename) :

    # Open the file

    file = h5py.File(filename, 'r')

    # Read attributes

    data = dict(zip(file.attrs.keys(),file.attrs.values()))

    # Read datasets

    for k in file.keys() :
        data[k] = file[k][...]

    # Convert items to complex

    complex_dtype = np.dtype([('re', '<f8'), ('im', '<f8')])

    for k in data.keys() :
        if(data[k].dtype == complex_dtype) :
            data[k] = data[k].astype(complex)

    # Return the data

    return data

#
