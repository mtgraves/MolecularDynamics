# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# This (will be) a class which holds multiple lattice structures.
# Right now it is a module of functions for getting an FCC lattice
# with N constituent atoms.
#
# Author:       Max Graves
# Date:         27-MAY-2013 (Memorial Day)
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import argparse, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

def plot3dim(X, Y, Z, color):
    ''' plot x,y,z components of atoms '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, s=60, c=color)
    ax.set_xlabel('$x\ (\AA)$')
    ax.set_ylabel('$y\ (\AA)$')
    ax.set_zlabel('$z\ (\AA)$')
    ax.set_title('FCC Lattice - Initial spatial configuration')
    plt.show()
    return 0

def plot2dim(X, Y, Z, color):
    ''' create 3 subplots of all three combinations of x,y,z plots '''
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

def buildUnitCell_FCC(a):
    """
    Build the unit cell for an FCC lattice.
    """
    uxpos = np.array([])
    uypos = np.array([])
    uzpos = np.array([])
    ucolor = np.array([])
    
    # define fcc vectors
    a_1 = 0.500*a*np.array([0,1,1])
    a_2 = 0.500*a*np.array([1,0,1])
    a_3 = 0.500*a*np.array([1,1,0])
    fccVecs = np.array([a_1, a_2, a_3])
    
    # define Tetrahedral basis, each should get a different color
    r_1 = 0.250*a*np.array([0,0,0])             # assign blue
    r_2 = 0.250*a*np.array([1,1,0])             # assign green
    r_3 = 0.250*a*np.array([1,0,1])             # assign red
    r_4 = 0.250*a*np.array([0,1,1])             # assign yellow
    tetBasis = np.array([r_1, r_2, r_3, r_4])   #[b, g, r, y]
    
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

    return uxpos, uypos, uzpos, ucolor

def tileEntireLatt(N, uzpos, uypos, uxpos, ucolor, a):
    """
    Tile out the lattice based on the unit vectors.
    """
    if N ==16:
        Nx, Ny, Nz = 1, 1, 1
    elif N == 128:
        Nx, Ny, Nz = 2, 2, 2

    count = 0
    xpos, ypos, zpos  = np.array([]), np.array([]), np.array([])
    color = np.array([])
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                for _ in range(uypos.size):
                    count += 1
                    if count == N+1:
                        print 'less atoms than desirable in cell'
                        return xpos, ypos, zpos, color
                    else:
                        xpos = np.append(xpos, uxpos[_]+i*a)
                        ypos = np.append(ypos, uypos[_]+j*a)
                        zpos = np.append(zpos, uzpos[_]+k*a)
                        color = np.append(color, ucolor[_])
 
    return xpos, ypos, zpos, color

def parseCMD():
    """ parse command line options """
    parser = argparse.ArgumentParser(description='MD for Argon')
    parser.add_argument("-N", "--numAtoms", type=int, default=128,
            help="number of Atoms to simulate (require Nmodn^3 = 0)")
    parser.add_argument("-a", "--lattSpacing", type=float, default=10,
            help="the Lattice constant in Angstroms")
    return parser.parse_args()
