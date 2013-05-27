# ===================================================================
# MD code for N non-interacting classical particles in 3 dimensions
#   using the Verlet algorithm for integration.
#
# Author:       Max Graves
# Date:         26-MAY-2013
# ===================================================================

import numpy as np
import pylab as pl

# ===================================================================
# harmonic oscillator potential and gradV
# ===================================================================
def Vharm(x,m):
    """ returns potential for harmonic oscillator """
    omega = 1.0
    return 0.5*m*omega**2*x**2

def gradVharm(x,m):
    """ returns gradient of potential for harmonic oscillator """
    omega = 1.0
    return m*omega*x

# ===================================================================
# main function
# ===================================================================
def main():

    N = 3       # number of particles
    D = 3       # number of dimensions
    m = 1.0     # mass of particles
    dt = 0.04   # time step
    
    r = np.array([[2.0, 1.5, 1.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]])
    v = np.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]])
    a = np.array(-gradVharm(r,m)/m)
    
    # array for time
    t = np.arange(dt,60,dt)
    numTimeSteps = t.size-1
    
    # running arrays to keep position, velocity, acceleration
    positions = np.array([r])
    velocitys = np.array([v])
    acceleras = np.array([a])

    # numerically integrate using the Verlet algorithm
    for i in range(numTimeSteps):
        v = v + 0.5*a*dt        # update half-step velocity
        r = r + dt*v            # update position
        a = -gradVharm(r,m)/m   # update acceleration
        v = v + 0.5*dt*a        # update velocity
        positions = np.append(positions, r)
        velocitys = np.append(velocitys, v)
        acceleras = np.append(acceleras, a)

    # reshape arrays for easy plotting of x,y,z coordinates
    pos = np.reshape(positions, ((numTimeSteps+1),D*N))
    vel = np.reshape(velocitys, ((numTimeSteps+1),D*N))
    acc = np.reshape(acceleras, ((numTimeSteps+1),D*N))

    pnew, vnew, anew  = np.array([]), np.array([]), np.array([])
    for i in range(len(pos[0])):
        for j in range(len(pos)):
            pnew = np.append(pnew, pos[j][i])
            vnew = np.append(vnew, vel[j][i])
            anew = np.append(anew, acc[j][i])
    pnew = np.reshape(pnew, (N, D, numTimeSteps+1))
    vnew = np.reshape(vnew, (N, D, numTimeSteps+1))
    anew = np.reshape(anew, (N, D, numTimeSteps+1))

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
