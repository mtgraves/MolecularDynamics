# =============================================================================
# Integration Methods.
#
# Author:               Max Graves
# Last Revision:        7-OCT-2012
# =============================================================================

import pylab as pl

# =============================================================================
class NewtonCotes:
    """
    Base class for all numerical integration methods which
    fall in the Newton-Cotes category.  Default number of panels
    is set to 1000, but can be over-ridden by simply passing 
    another N value from your program.
    """
    def __init__(self, f, a, b, N=1000):
        """
        Takes variables as described in comments below.
        This sets up variables needed in Newton-Cotes
        integration schemes.
        """
        self.dx = (b-a)/N   # spatial discretization
        self.I = 0.0        # value of integral
        self.f = f          # function
        self.a = a          # lower integration limit
        self.b = b          # upper integration limit
        self.N = N          # number of panels

        # instantaneous values of integral (for plotting)
        self.plot = pl.array([])

    def DiscreteSum(self,lower,upper):
        """
        Perform the rectangular integration, starting
        and ending at the values passed to it (0-->N for rect,
        1-->N for trap, etc...).
        """
        for i in range(int(lower),int(upper)):
            self.I += self.f(self.a+i*self.dx)*self.dx
            self.plot = pl.append(self.plot, self.I)
        return self.I, self.plot

# =============================================================================
class Rectangular(NewtonCotes):
    """
    Perform Rectangular method of numerical integration.
    This is a Newton-Cotes method.
    Inherits from base class.
    """
    def rect(self):
        """
        Rectangular numerical quadrature. Returns the 
        value of the definite integral as well as an array for
        plotting of the instantaneous integral values over the 
        integration variable.
        """
        self.I, self.plot = self.DiscreteSum(0,self.N)
        return self.I, self.plot

# =============================================================================
class Trapezoidal(NewtonCotes):
    """
    Perform Trapezoidal method of numerical integration.
    This is a Newton-Cotes method.
    Inherits from base class.
    """
    def trap(self):
        """
        Trapezoidal numerical quadrature. Returns the 
        value of the definite integral as well as an array for
        plotting of the instantaneous integral values over the 
        integration variable.
        """
        self.I += (self.dx/2.0)*self.f(self.a)
        self.plot = pl.append(self.plot,self.I)
        self.I,self.plot = self.DiscreteSum(1,self.N)
        self.I += (self.dx/2.0)*self.f(self.b)
        return self.I, self.plot

# =============================================================================
class Romberg(NewtonCotes):
    """
    Perform Romberg method of numerical integration.
    This recursively calculates better approximations to the
    value of the integral, and will calculate the terms of the 
    trace of this lower triangular matrix up to 6 decimal places
    of precision.  Inherits from base class.
    """
    def romb(self):
        """
        Romberg numerical quadrature to 6 digits of agreement
        between subsequent diagonal elements.
        """
        R = pl.zeros([1,1])
        # check convergence twice (to avoid flukes)
        Rdiff1 = 100.0
        Rdiff2 = 100.0
        # zeroeth order term
        constant = 0.50*(self.b-self.a)
        R[0,0] = constant*(self.f(self.a)+self.f(self.b))
        count = 1
        # higher order terms
        while (Rdiff2 > 0.000001) or (Rdiff1 > 0.000001) :
            
            # increase size of matrix by 1 row,column
            R = pl.append(R, pl.zeros([1,count]), 0)
            R = pl.append(R, pl.zeros([count+1,1]), 1)

            # build next row, and only that row
            for n in range(count,count+1):
                dx = (self.b-self.a)/(2.0**n)
                for m in range(count+1):
                    if m<=n:
                        if m==0:
                            R[n,m] += 0.5*R[n-1,m]
                            for j in range(1,int(2**(n-1)+1)):
                                R[n,m] += dx*self.f(self.a+(2*j-1)*dx)
                        else:
                            term1 = R[n,m-1]
                            term2 = R[n-1,m-1]
                            R[n,m] = ((4**m)*term1-term2)/(4**m-1)
            
            count += 1
            
            Rdiff1 = abs(R[-1,-1] - R[-2,-2])
            if count > 3:
                Rdiff2 = abs(R[-2,-2] - R[-3,-3])
        return R[-1,-1]

# =============================================================================
