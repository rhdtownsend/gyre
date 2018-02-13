#!/usr/bin/env python
#
# Convert an HDF-format file to a TXT-format file

import gyre
import sys
import numpy as np
import fortranformat as ff

# Arguments

in_file = sys.argv[1]
out_file = sys.argv[2]

# Formatted writers

ind_form = ff.FortranRecordWriter('(I25)')
name_form = ff.FortranRecordWriter('(A25)')

int_form = ff.FortranRecordWriter('(I25)')
float_form = ff.FortranRecordWriter('(E25.16E3)')

# Read the GYRE data

d, r = gyre.read_output(in_file)

# Write the output data

with open(out_file, 'w') as f:

    # Label

    f.write(d['label'].flatten()[0].decode('UTF-8')+'\n')

    # Attributes

    i = 0

    ind_str = ''
    name_str = ''
    val_str = ''

    for name in d.dtype.names:

        if name == 'label':
            continue

        dtype = d[name].dtype

        if dtype == 'int32':

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write([name])
            val_str += int_form.write([d[name].flatten()[0]])

        elif dtype == 'float32' or dtype == 'float64':

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write([name])
            val_str += float_form.write([d[name].flatten()[0]])

        elif dtype == 'complex64' or dtype == 'complex128':

            i += 1

            ind_str += int_form.write([i])
            name_str += name_form.write([name])
            val_str += float_form.write([d[name].flatten()[0].real])

            i += 1

            ind_str += int_form.write([i])
            name_str += name_form([name])
            val_str += float_form.write([d[name].flatten()[0].imag])

        else:

            raise Exception('Cant deal with field: {:s}'.format(name))

    f.write(ind_str+'\n')
    f.write(name_str+'\n')
    f.write(val_str+'\n')

    # Datasets
                        
    i = 0
    n = 0

    ind_str = ''
    name_str = ''
    val_str = ''

    for name in r.dtype.names:

        if name == 'i':
            continue

        if n == 0:
            n = len(r[name])
            val_str = np.empty([n], dtype=object)
            for j in range(n):
                val_str[j] = ''
            
        dtype = r[name].dtype

        if dtype == 'int32':

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write([name])
            for j in range(n):
                val_str[j] += int_form.write([r[name][j]])

        elif dtype == 'float32' or dtype == 'float64':

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write([name])
            for j in range(n):
                val_str[j] += float_form.write([r[name][j]])

        elif dtype == 'complex64' or dtype == 'complex128':

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write(['Re({:s})'.format(name)])
            for j in range(n):
                val_str[j] += float_form.write([r[name][j].real])

            i += 1

            ind_str += ind_form.write([i])
            name_str += name_form.write(['Im({:s})'.format(name)])
            for j in range(n):
                val_str[j] += float_form.write([r[name][j].imag])

        else:

            raise Exception('Cant deal with field: {:s}'.format(name))

    f.write(ind_str+'\n')
    f.write(name_str+'\n')
    for j in range(n):
        f.write(val_str[j]+'\n')

        


        
