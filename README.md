MolecularDynamics
=================

Classical MD under development, soon to be extended to PIMD.

Uses the Verlet Velocity algorithm to integrate the equations to motion.
This is supposed to be good to fourth order in the time step.

Done
=======

- Lennard-Jones 6-12 potential and gradient implemented
    - currently counting all two-body interactions in the slowest way possible

To Do
======

- Periodic Boundary conditions
- ability to start atoms in lattice configuration (fcc to begin with)
