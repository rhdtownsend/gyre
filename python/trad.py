import numpy as np
import numpy.polynomial.chebyshev as cheby

class Trad :

    """Traditional-approximation data"""

    def __init__ (self, m, k, cb_pos, cb_neg, cb_ctr) :

        assert isinstance(cb_pos, cheby.Chebyshev)
        assert isinstance(cb_ctr, cheby.Chebyshev)
        assert isinstance(cb_neg, cheby.Chebyshev)

        self.cb_pos = cb_pos
        self.cb_neg = cb_neg
        self.cb_ctr = cb_ctr

        self.m = m
        self.k = k

    def l_e (self, nu) :

        """Effective harmonic degree"""

        l_e = 0.5*(-1. + np.sqrt(1. + 4.*self.lamb(nu)))

        return l_e

    def lamb (self, nu) :

        """Laplace tidal equation eigenvalue"""

        if isinstance(nu, float) :
            nu_r = nu
        elif isinstance(nu, complex) :
            nu_r = nu.real
        else :
            raise Exception('nu must be float or complex')

        if nu_r < -1. :
            lamb = self.cb_neg(1./nu)*self.lamb_asymp(nu)
        elif nu_r > 1. :
            lamb = self.cb_pos(1./nu)*self.lamb_asymp(nu)
        else :
            l = abs(self.m) + self.k
            lamb = self.cb_ctr(nu)*l*(l+1)

        return lamb

    def lamb_asymp (self, nu) :

        """Laplace tidal equation eigenvalue in asymptotic limit of large nu"""

        if isinstance(nu, float) :
            nu_r = nu
        elif isinstance(nu, complex) : 
            nu_r = nu.real
        else :
            raise Exception('nu must be float or complex')
        
        if self.m*nu_r >= 0. :

            assert self.k >= 0

            if self.k > 0 : 
                s = self.k - 1
                lamb_asymp = self.m*nu + self.m**2 + 0.5*nu**2*(2*s + 1)**2*(1. + np.sqrt(1. + 4.*(self.m*nu + self.m**2)/(nu**2*(2*s + 1)**2)))
            else :
                lamb_asymp = self.m**2*(2.*self.m*nu)/(2.*self.m*nu - 1.)

        else :

            if self.k >= -1 : 
                s = self.k + 1
                lamb_asymp = self.m*nu + self.m**2 + 0.5*nu**2*(2*s + 1)**2*(1. + np.sqrt(1. + 4.*(self.m*nu + self.m**2)/(nu**2*(2*s + 1)**2)))
            else :
                s = -self.k -1
                lamb_asymp = (self.m*nu - self.m**2)**2/(nu**2*(2*s+1)**2)

        return lamb_asymp
