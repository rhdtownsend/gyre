#!/usr/bin/env python
#
# Build reference files

import numpy as np

# Build a reference file using analytical frequencies (from [Pek1938])

def build_ref_0 (out_file, gamma) :

    # Write the output file

    f = open(out_file, 'w')

    f.write("\n\n\n\n")

    f.write("1 2 3 4 5 6 7\n")
    f.write("l n_pg n_p n_g Re(omega) Im(omega) E_norm\n")
    
    for l in range(0, 4) :

        for n in range(1, 11) :

            k = 2*(n-1)

            D = 0.5*(-4. + 0.5*gamma*(k*(k+5+2*l) + 6 + 4*l))

            beta = D + np.sqrt(D**2 + l*(l+1))

            omega = np.sqrt(beta)

            f.write("{:d} {:d} 0 0 {:12.9f} 0. 0.\n".format(l, n, omega))

    f.close()
            

# Build a reference file using tabulated frequency data

def build_ref (in_file, out_file, nu_ref, l_0) :

    # Read the input file

    d = np.loadtxt(in_file)

    # Write the output file

    f = open(out_file, 'w')

    f.write("\n\n\n\n")

    f.write("1 2 3 4 5 6 7\n")
    f.write("l n_pg n_p n_g Re(omega) Im(omega) E_norm\n")
    
    for j in range(1, d.shape[1]) :

        l = j - 1 + l_0

        for i in range(0, d.shape[0]) :

            n_pg = int(np.rint(d[i,0]))
            nu = d[i,j]

            f.write("{:d} {:d} 0 0 {:12.9f} 0. 0.\n".format(l, n_pg, nu/nu_ref))

    f.close()

#

gamma = 5./3.

build_ref_0('0.0/ref/summary.txt', gamma)

nu_ref = 99.855377

build_ref('data/poly-freqs-1.5.dat', '1.5/ref/summary.txt', nu_ref, 0)
build_ref('data/poly-freqs-3.0.dat', '3.0/ref/summary.txt', nu_ref, 0)
build_ref('data/poly-freqs-3.0-g.dat', '3.0-g/ref/summary.txt', nu_ref, 1)
build_ref('data/poly-freqs-4.0.dat', '4.0/ref/summary.txt', nu_ref, 0)
