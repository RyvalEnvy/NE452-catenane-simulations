from openmm import app
import openmm as mm
from openmm import unit
import sys
import h5py
import numpy as np
import os

def runSim(T = 50, gamma = 1/50, len = 10, dt = 0.1, skipSteps = 1, f_pdb = "pdb", 
           f_ff = "forcefield"):
    """
    Function to run the openMM simulation (NVT) with given params and 
    pdb+forcefields 

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
    """
    # set up strings
    # THE FOLLOWING FILE NAMES ARE PLACEHOLDERS: CHANGE WHEN POSSIBLE
    path = os.getcwd()
    pdbPath = os.path.join(path, "pdb_files/{}.pdb".format(f_pdb))
    forcefieldPath = os.path.join(path, "pdb_files/{}.xml".format(f_ff))

    beta = 1./T
    steps = int(len/dt)
    indSteps = int(len/(dt*skipSteps))

    pdb = app.PDBFile(pdbPath)

    topo = pdb.topology

    particles = []
    for atom in topo.atoms():
        particles.append(atom)

    # rough shitty code below NOT VERY ROBUST
    carbons = particles[0:20]
    hy1 = particles[20:40]
    hy2 = particles[40:60]

    for i, C in enumerate(carbons):
        # print("Carbon: {}, {}".format(C, Cs[i-1]))
        # print("Hydrogen: {}, {}".format(H1s[i], H2s[i]))
        topo.addBond(C, carbons[i-1])
        topo.addBond(C, hy1[i])
        topo.addBond(C, hy2[i])

    forcefield = app.ForceField(forcefieldPath)
    unmatched_residues = forcefield.getUnmatchedResidues(topo)
    print("unmatched residues\n", unmatched_residues)
    nonbonded = app.NoCutoff
         
    system = forcefield.createSystem(topo, nonbondedMethod = nonbonded,
                                    nonbondedCutoff = 1e3*unit.nanometer,
                                    constraints = None)

    integrator = mm.LangevinIntegrator(T*unit.kelvin, gamma, dt*unit.femtoseconds)

    platform = mm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.computeVirtualSites()
    	
    simulation.context.setVelocitiesToTemperature(T*unit.kelvin)

    state = simulation.context.getState(getForces = True, getEnergy = True,
                                        getPositions = True)
    
    # generate h5 to store info about positions of nodes in sim
    # outH5 = h5py.File("cattraj_{T}K.h5".format(T=T))
    # h5coords=outH5.create_dataset("coordinates",(int(steps/skipSteps)+1,
    #                                 nAtoms,3), dtype = "f")
    # h5kEnergy=outH5.create_dataset("kineticEnergy",(int(steps/skipSteps)+1), 
    #                                 dtype = "f")
    # h5pEnergy=outH5.create_dataset("potentialEnergy",(int(steps/skipSteps)+1), 
    #                                 dtype = "f")
    # h5coords.attrs['units'] = "nanometers"
    # h5coords.attrs['frame'] = "lab-fixed"
    # h5coords.attrs['temperature'] = T
    # h5coords.attrs['steps'] = steps
    # h5coords.attrs['skipSteps'] = skipSteps
    # h5coords.attrs['dt'] = dt/unit.femtosecond
    # h5coords.attrs['gamma'] = gamma

    for step in range(steps):
        simulation.step(skipSteps)
        state = simulation.context.getState(getEnergy=True)
        KE = state.getKineticEnergy()/unit.kilojoules_per_mole
        PE = state.getPotentialEnergy()/unit.kilojoules_per_mole
        print(KE)
        print(PE)

if __name__ == "__main__":
     runSim(T = 50, gamma = 1/50, len = 1, dt = 0.1, skipSteps = 1, f_pdb = "small_ring_full", 
           f_ff = "forcefield_sm")