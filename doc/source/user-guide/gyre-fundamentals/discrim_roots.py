#!/usr/bin/env python
#
# Discriminant function roots calculator

# Imports

import numpy as np
import numpy.linalg as la
import scipy.optimize as op

from colors import *

# Calculation parameters

N = 50

sigma_min = 0
sigma_max = 5.5
delta_sigma = 0.17

# System matrix

def S (sigma):

    tau = np.pi/(N-1)

    S = np.zeros([N, N])

    S[0,0] = 1.

    for i in range(1, N-1):
        S[i,i] = tau**2*sigma**2 - 2.
        S[i,i-1] = 1.
        S[i,i+1] = 1.

    S[N-1,N-1] = 1.

    return S

# Discriminant function

def discrim (sigma):

    return la.det(S(sigma))/((-1)**(N-2)*(N-1))

# Evaluate discriminant data

sigma = np.concatenate([np.arange(0, sigma_min, -delta_sigma)[1:], np.arange(0, sigma_max, delta_sigma)])
n_sigma = len(sigma)
D = np.empty(n_sigma)

for i in range(n_sigma):
    D[i] = discrim(sigma[i])

# Bracket and find roots

i_bracket = np.where(D[1:]*D[:-1] <= 0.)[0]
n_bracket = len(i_bracket)

sigma_root = np.empty(n_bracket)

for j in range(n_bracket):
    sigma_root[j] = op.brentq(discrim, sigma[i_bracket[j]], sigma[i_bracket[j]+1])

# Write out roots in CSV format

with open('eigenfreqs.csv', 'w') as f:

    f.write('"n", "numerical", "analytic"\n')

    for j in range(n_bracket):
        f.write("{:d}, {:8.6f}, {:8.6f}\n".format(j+1, sigma_root[j], float(j+1)))
