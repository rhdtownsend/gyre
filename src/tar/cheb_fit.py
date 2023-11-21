# Python class implementing GYRE's cheby_fit_t

import numpy as np
import h5py

class ChebFit:

    def __init__ (self, x_a, x_b, f):

        self.x_a = x_a
        self.x_b = x_b

        self.f = f

    @classmethod
    def load (cls, filename, group='/'):

        with h5py.File(filename, 'r') as f:

            x_a = f[group].attrs['x_a']
            x_b = f[group].attrs['x_b']

            f = f[group]['f'][...]

        return cls(x_a, x_b, f)

    def eval (self, x):

        if x == self.x_a:
            u = 1.
        elif x == self.x_b:
            u = -1.
        else:
            u = (2.*x - (self.x_a + self.x_b))/(self.x_a - self.x_b)

        s_n = 0.
        s_d = 0.

        n = len(self.f) - 1

        for j in range (0, n+1):

            u_j = np.cos(j*np.pi/n)

            if u == u_j:
                return self.f[j]

            if j == 0 or j == n:
                w = 0.5*(-1.)**j
            else:
                w = (-1.)**j

            s_n += w*self.f[j]/(u - u_j)
            s_d += w/(u - u_j)

        return s_n/s_d

    def coeffs (self):

        n = len(self.f) - 1

        c = np.empty([n+1])
   
        for k in range(0, n+1):

            c[k] = 0.

            for j in range(0, n+1):

                if j == 0:
                    v = 1.
                    w = 0.5
                elif j == n:
                    v = (-1.)**k
                    w = 0.5
                else:
                    v = np.cos(j*k*np.pi/n)
                    w = 1.

                c[k] += w*v*self.f[j]

            if k == 0 or k == n:
                c[k] /= n
            else:
                c[k] /= 0.5*n
        
        return c
    
