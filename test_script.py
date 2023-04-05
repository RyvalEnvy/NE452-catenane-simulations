from openmm import app
import openmm as mm
from openmm import unit
import sys
import h5py
import numpy as np
import os

def runSim(T = 50, gamma = 1/50, length = 10, dt = 0.1, skipSteps = 1, f_pdb = "pdb", 
           f_ff = "forcefield", srSize = 20):
    """
    Function to run the openMM simulation (NVT) with given params and 
    pdb+forcefields

    Saves a xyz file of a compressed 3d array. See the function reload() to get 
    an idea of how to load and manipulate the output array (and format the file 
    name). The final array should have dims defined by 
    [n_timesteps, 3 (xyz), n_particles]

    parameters:
    T: int
        Temperature
    gamma: float
        Friction parameters for Langevin
    len: float
        Total length of the simulation in ps
    dt: float
        Size of the timesteps in ps
    f_pdb: str
        Name of the pdb file (excluding the file extension, which should 
        simply be .pdb)
    f_ff: str
        Name of the forcefield file (excluding the file extension, which should 
        simply be .xml)
    srSize: int
        Number of monomer units in the small ring
    """
    # set up strings
    # THE FOLLOWING FILE NAMES ARE PLACEHOLDERS: CHANGE WHEN POSSIBLE
    path = os.getcwd()
    pdbPath = os.path.join(path, "pdb_files/{}.pdb".format(f_pdb))
    forcefieldPath = os.path.join(path, "pdb_files/{}.xml".format(f_ff))

    # set up params
    kB = 1.38e-23
    beta = 1./(kB * T) # needs kb
    steps = int(length*1000/dt)
    indSteps = int(length*1000/(dt*skipSteps))

    # load external files for topo and PES
    pdb = app.PDBFile(pdbPath)
    topo = pdb.topology

    # identify atoms in topology
    particles = []
    # future update maybe: collect residues to identify number of atoms/chain
    # residues = []
    for atom in topo.atoms():
        particles.append(atom)
        # residues.append(atom.residue)

    # add bonds between the individual rings
    c1 = particles[0:srSize]
    c2 = particles[srSize:]
    numParticles = len(particles)
    # hy1 = particles[20:40]
    # hy2 = particles[40:60]

    # create arrays to store xyz coords in
    # format is timestep:xyz:particle#
    traj = np.zeros((indSteps, 3, numParticles))

    for i, C in enumerate(c1):
        topo.addBond(C, c1[i-1])
        # topo.addBond(C, hy1[i])
        # topo.addBond(C, hy2[i])
    
    for i, C in enumerate(c2):
        topo.addBond(C, c2[i-1])
        # topo.addBond(C, hy1[i])
        # topo.addBond(C, hy2[i])

    forcefield = app.ForceField(forcefieldPath)
    unmatched_residues = forcefield.getUnmatchedResidues(topo)
    print("unmatched residues\n", unmatched_residues)
    nonbonded = app.NoCutoff
         
    system = forcefield.createSystem(topo, nonbondedMethod = nonbonded,
                                    nonbondedCutoff = 1e3*unit.nanometer,
                                    constraints = None)

    integrator = mm.LangevinIntegrator(T*unit.kelvin, gamma, 
                                       dt*unit.femtoseconds)

    platform = mm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.computeVirtualSites()
    	
    simulation.context.setVelocitiesToTemperature(T*unit.kelvin)

    # equilibration step
    # simulation.step(1000)

    for step in range(steps):
        simulation.step(skipSteps)
        # extract positions
        current = simulation.context.getState(getPositions = True)
        positions = np.array(current.getPositions()/unit.nanometer)
        traj[step, :, :] = positions.transpose()

    # save the trajectory as an xyz file
    traj_r = traj.reshape(traj.shape[0], -1)
    print(traj_r.shape)
    np.savetxt("{}fs_{}ps_{}K_{}ss.xyz".format(dt, length, T, skipSteps), traj_r)

def reload(T = 50, gamma = 1/50, length = 10, dt = 0.1, skipSteps = 1, 
           f_pdb = "pdb", f_ff = "forcefield", srSize = 20):
    fname = "{}fs_{}ps_{}K_{}ss.xyz".format(dt, length, T, skipSteps)
    totSteps = int(length*1000/(skipSteps * dt))
    print("totSteps", totSteps)

    a = np.loadtxt(fname)
    nParticles = a.shape[-1]//3
    a = a.reshape(totSteps, 3, nParticles)
    print(a.shape)

    return a

if __name__ == "__main__":
    runSim(T = 50, gamma = 1/50, length = 0.01, dt = 0.1, skipSteps = 1, 
           f_pdb = "combo_ring", f_ff = "combo_ff")
    # reload(T = 50, gamma = 1/50, length = 0.01, dt = 0.1, skipSteps = 1, 
    #        f_pdb = "combo_ring", f_ff = "combo_ff")