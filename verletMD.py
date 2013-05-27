# ===================================================================
# MD code for N Argon atoms in 3-D, using a Lennard-Jones potential.
# Right now, this implements the velocity-Verlet algorithm.
#
#   V(r) = 4e( (s/r)^12 - (s/r)^6 )
#
# Argon specs:
#       m   = 6.69*10^-26 kg
#       e   = 1.654*10^-21 J
#       s   = 3.405*10^-10 m
#       r_m = 2^(1/6)*s
#
#       one unit of time corresponds to tau = 2.17 ps
#
# Author:       Max Graves
# Date:         26-MAY-2013
# ===================================================================

import numpy as np
import pylab as pl
import argparse, sys, fccLatt, random

def LennardJones(r, i):
    """ 
    Compute the potential energy for the i^th particle in 
    a system using the Lennard-Jones 6-12 potential.  Also compute
    the gradient of the potential.  Returns both.
    """
    s = 1.0             # sigma (see header line)
    e = 1.0             # epsilon (see header line)
    V = 0.0
    gradV = np.array([0,0,0])
    # r[i][n] gives access to x=0,y=1,z=2 coordinate of i^th atom
    for a in range(len(r)):
        if a == i:
            pass
        else:
            R = 0.0
            for n in range(len(r[0])):
                R += (r[a][n] - r[i][n])**2
            R = np.sqrt(R)
            V += 4.0*e*((1.0*s/R)**12 - (1.0*s/R)**6)
            gradV = np.add(gradV, 24.0*e*(r[a] - r[i])*( 2.0*s**12/(R**14) -
                s**6/(R**8)))

    return V, gradV

def parseCmd():
    """ parse command line options """
    parser = argparse.ArgumentParser(description='MD for Argon')
    parser.add_argument("-N", "--numAtoms", type=int, default=128,
            help="Number of atoms to simulate")
    parser.add_argument("-a", "--lattSpacing", type=float, default=10,
            help="Lattice spacing in Angstroms")
    parser.add_argument("-T", "--temp", type=float, default=70,
            help="Temperature in Kelvin")
    parser.add_argument("-D", "--dimension", type=int, default=3,
            help="Number of dimensions")
    return parser.parse_args()

# =============================================================================
def main():

    args = parseCmd()

    N = args.numAtoms       # number of particles
    D = args.dimension      # number of dimensions
    a = args.lattSpacing    # lattice constant (Angstroms)
    T = args.temp
    m = 1.0                 # mass (see header line)
    dt = 0.04               # time step
    k = 1.3806488*10**(-23) # Boltzmann constant
  
    # set initial positions as FCC lattice
    uxpos, uypos, uzpos, ucolor = fccLatt.buildUnitCell_FCC(a)

    xpos, ypos, zpos, color = fccLatt.tileEntireLatt(N,
            uzpos, uypos, uxpos, ucolor, a)

    r = np.array([])
    for i in range(len(xpos)):
        r = np.append(r, [xpos[i], ypos[i], zpos[i]])
    r = np.reshape(r, (N,D))
    
    # set initial velocities randomly from Boltzmann distribution
    v = np.array([])
    for i in range(N*D):
        ran = random.gauss(0,1)
        v = np.append(v, 248.822*ran)   # v = sqrt(kT/m)*ran
    v = np.reshape(v, (N,D))

    # compute acceleration due to LJ potential
    a = np.array([])
    for i in range(N):
        a = np.append(a, -(1.0/m)*LennardJones(r,i)[1])
    a = np.reshape(a,(N,D))

    # array for time
    t = np.arange(dt,60,dt)
    numTimeSteps = t.size-1
    
    # running arrays to keep position, velocity, acceleration
    positions = np.array([r])

    # numerically integrate using the Verlet algorithm
    for i in range(numTimeSteps):
            v = v + 0.5*a*dt        # update half-step velocity
            r = r + dt*v            # update position
            for j in range(N):
                a[j] = -(1.0/m)*LennardJones(r,j)[1]   # update acceleration
            v = v + 0.5*dt*a        # update velocity
            positions = np.append(positions, r)
    
    # reshape arrays for easy plotting of x,y,z coordinates
    pos = np.reshape(positions, ((numTimeSteps+1),D*N))
    pnew  = np.array([])
    for i in range(len(pos[0])):
        for j in range(len(pos)):
            pnew = np.append(pnew, pos[j][i])
    pnew = np.reshape(pnew, (N, D, numTimeSteps+1))

    # Plot
    pl.plot(t, pnew[0][0], color='k', label='x component')
    pl.plot(t, pnew[0][1], color='Lime', label='y component')
    pl.plot(t, pnew[0][2], color='Violet', label='z component')
    pl.xlabel('time '+r'$[s]$')
    pl.ylabel('position '+r'$[m]$')
    pl.legend()
    pl.grid(True)
    pl.show()

# ===================================================================
if __name__=="__main__":
    main()
