# =============================================================================
# Class for defining lattices (INCOMPLETE!)
#
# Author:               Max Graves
# Last Revision:        27-MAY-2013 (Memorial day)
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
class Cubic:
    """
    Base class for all cubic lattice structures.
    """
    def __init__(self, f, a, Nx, Ny, Nz):
        """
        Takes variables as described in comments below.
        Sets up variables for lattice as described in
        comments below.
        """
        self.a = a          # lattice constant
        self.Nx = Nx        # Number of tiles in x
        self.Ny = Ny        # Number of tiles in y
        self.Nz = Nz        # Number of tiles in z

        # define primitive lattice vectors
        self.a_1 = np.array([])
        self.a_2 = np.array([])
        self.a_3 = np.array([])
        self.primVecs = np.array(
                [self.a_1, self.a_2, self.a_3])

        # define basis vectors
        self.r_1 = np.array([])
        self.r_2 = np.array([])
        self.r_3 = np.array([])
        self.r_4 = np.array([])
        self.BasisVecs = np.array(
                [self.r_1, self.r_2, self.r_3, self.r_4])

        # arrays for unit cell
        self.uxpos = np.array([])
        self.uypos = np.array([])
        self.uzpos = np.array([])
        self.ucolor = np.array([])

        # arrays for entire lattice
        self.xpos = np.array([])
        self.ypos = np.array([])
        self.zpos = np.array([])
        self.color = np.array([])

    def buildUnit(self):
        """
        Tile out the unit cell.  This should
        be very general and should work no matter
        what the lattice as long as it is monatomic.
        """
        for m in range(tetBasis.size/3):
        uxpos = np.append(uxpos, tetBasis[m][0])
        uypos = np.append(uypos, tetBasis[m][1])
        uzpos = np.append(uzpos, tetBasis[m][2])
        if m == 0:
            ucolor = np.append(ucolor, 'b')
        elif m == 1:
            ucolor = np.append(ucolor, 'g')
        elif m == 2:
            ucolor = np.append(ucolor, 'r')
        elif m == 3:
            ucolor = np.append(ucolor, 'y')
        for n in range(fccVecs.size/3):
            uxpos = np.append(uxpos, tetBasis[m][0]+fccVecs[n][0])
            uypos = np.append(uypos, tetBasis[m][1]+fccVecs[n][1])
            uzpos = np.append(uzpos, tetBasis[m][2]+fccVecs[n][2])
            if m == 0:
                ucolor = np.append(ucolor, 'b')
            elif m == 1:
                ucolor = np.append(ucolor, 'g')
            elif m == 2:
                ucolor = np.append(ucolor, 'r')
            elif m == 3:
                ucolor = np.append(ucolor, 'y')


    
    def plot3dim(self, X, Y, Z, color):
        """
        plot x,y,z components of atoms
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X, Y, Z, s=60, c=color)
        ax.set_xlabel('$x\ (\AA)$')
        ax.set_ylabel('$y\ (\AA)$')
        ax.set_zlabel('$z\ (\AA)$')
        ax.set_title('Pyrochlore Lattice Unit cell')
        plt.show()
        return 0

    def plot2dim(self, X, Y, Z, color):
        """
        Create 3 subplots of all three combinations of x,y,z plots
        """
        fig = plt.figure()
        gs = gridspec.GridSpec(4,4)

        xyplot = plt.subplot(gs[0:2,:2])
        xyplot.set_title(r'$x,y\ positions$')
        xyplot.scatter(X,Y, c=color)
        xyplot.set_ylabel(r'$y\ (\AA)$')
        xyplot.set_xlabel(r'$x\ (\AA)$')

        xzplot = plt.subplot(gs[2:4,:2])
        xzplot.set_title(r'$x,z\ positions$')
        xzplot.scatter(X,Z,c=color)
        xzplot.set_ylabel(r'$z\ (\AA)$')
        xzplot.set_xlabel(r'$x\ (\AA)$')

        yzplot = plt.subplot(gs[1:3,2:4])
        yzplot.set_title(r'$y,z\ positions$')
        yzplot.scatter(Y,Z, c=color)
        yzplot.set_ylabel(r'$z\ (\AA)$')
        yzplot.set_xlabel(r'$y\ (\AA)$')

        fig.tight_layout()
        plt.show()
        return 0

# =============================================================================
class FCC(Cubic):
    """
    Face-Centered-Cubic lattice.
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
class Pyrochlore(Cubic):
    """
    Pyrochlore lattice.
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
#class Romberg(NewtonCotes):
#    """
#    Perform Romberg method of numerical integration.
#    This recursively calculates better approximations to the
#    value of the integral, and will calculate the terms of the 
#    trace of this lower triangular matrix up to 6 decimal places
#    of precision.  Inherits from base class.
#    """
#    def romb(self):
#        """
#        Romberg numerical quadrature to 6 digits of agreement
#        between subsequent diagonal elements.
#        """
#        R = pl.zeros([1,1])
#        # check convergence twice (to avoid flukes)
#        Rdiff1 = 100.0
#        Rdiff2 = 100.0
#        # zeroeth order term
#        constant = 0.50*(self.b-self.a)
#        R[0,0] = constant*(self.f(self.a)+self.f(self.b))
#        count = 1
#        # higher order terms
#        while (Rdiff2 > 0.000001) or (Rdiff1 > 0.000001) :
#            
#            # increase size of matrix by 1 row,column
#            R = pl.append(R, pl.zeros([1,count]), 0)
#            R = pl.append(R, pl.zeros([count+1,1]), 1)
#
#            # build next row, and only that row
#            for n in range(count,count+1):
#                dx = (self.b-self.a)/(2.0**n)
#                for m in range(count+1):
#                    if m<=n:
#                        if m==0:
#                            R[n,m] += 0.5*R[n-1,m]
#                            for j in range(1,int(2**(n-1)+1)):
#                                R[n,m] += dx*self.f(self.a+(2*j-1)*dx)
#                        else:
#                            term1 = R[n,m-1]
#                            term2 = R[n-1,m-1]
#                            R[n,m] = ((4**m)*term1-term2)/(4**m-1)
#            
#            count += 1
#            
#            Rdiff1 = abs(R[-1,-1] - R[-2,-2])
#            if count > 3:
#                Rdiff2 = abs(R[-2,-2] - R[-3,-3])
#        return R[-1,-1]

# =============================================================================
