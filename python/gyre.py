# Module   : gyre
# Purpose  : process GYRE output
#
# Copyright 2013 Rich Townsend, Pieter Degroote
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

def read_output(filename):
    """
    Read a GYRE HDF5-format output file.
    
    Parameters
    ----------
    
    filename : string
        Input file
    
    Returns
    -------
    data: dictionary
        Dictionary containing the scalar information from the GYRE file. In the
        case of eigenvalue files, this dictionary is empty.
    local: record array, 
        Record array containing array information with the column names from
        the GYRE file. For eigenvalue files, this contains all the frequency
        information. For eigenfunction files, this contains all the
        eigenfunctions.
        
    Examples
    --------
    >>> data,local = gyre.read_output('eigvals.h5')
    
    You can access information by column:
    
    >>> local['omega']
    array([  2.01844773+0.j,   3.51781378+0.j,   5.02055334+0.j,
             6.38583592+0.j,   7.93991267+0.j,   9.50455841+0.j,
            11.01185876+0.j,  12.63161261+0.j,  14.26029392+0.j,
            15.84866171+0.j,  17.43948807+0.j,  19.04329558+0.j,
            20.66636443+0.j,  22.29938853+0.j,  23.91005215+0.j,
            25.52563311+0.j,  27.14974168+0.j])
    >>> local['n_g']
    array([ 3,  4,  5,  6,  7,  7,  9, 11, 11, 12, 13, 14, 15, 16, 17, 18, 19], dtype=int32)
    
    But also by row:
    
    >>> local[0]
    (1.0, 3, 0, (0.35264730283221724+0j), 2, (2.018447726710389+0j))        
    
    You can easily select for example all frequencies above a certain threshold:
    
    >>> selection = local[ (local['omega'].real > 10.) ]
    >>> selection['n_p']
    array([ 9, 11, 11, 12, 13, 14, 15, 16, 17, 18, 19], dtype=int32)
    
    If you just want to get information on the names of the columns, simply do:
    
    >>> local.dtype.names
    """

    # Open the file

    file = h5py.File(filename, 'r')

    # Read attributes

    data = dict(zip(file.attrs.keys(),file.attrs.values()))

    # Read datasets

    for k in file.keys() :
        data[k] = file[k][...]

    # Close the file

    file.close()

    complex_dtype = np.dtype([('re', '<f8'), ('im', '<f8')])
    
    # Convert the data into a record array and the global information in a dict
    
    local = []
    local_names = []
    
    for k in data.keys():
        
        # Convert items to complex
        
        if(data[k].dtype == complex_dtype) :
            data[k] = data[k]['re'] + 1j*data[k]['im']
        
        # Save array information to local record array, keep scalar information
        
        if not np.isscalar(data[k]):
            local.append(data.pop(k))
            local_names.append(k)
            
    local = np.rec.fromarrays(local,names=local_names)
    
    # Return the global and local information

    return data, local

#

if __name__=="__main__":
    import sys
    
    # Check if matplotlib is installed
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is not installed, so you have no access to the command line interface")
    
    # Read the info
    
    header,local = read_output(sys.argv[1])
    
    # If there is global information, print it out first
    if header:
        for key in header:
            print("{:10s} = {}".format(key,header[key]))
        print_resolution = len(local)/20
        print("Only showing 1 in every {} points in following array information".format(print_resolution))
        local = local[::print_resolution]
    
    # then print the local information
    
    print(plt.mlab.rec2txt(local))
