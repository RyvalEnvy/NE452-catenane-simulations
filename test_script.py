from openmm.app import app
from openmm import mm
from openmm import unit
import sys
import h5py
import numpy as np
import os

def runSim(T = 50, gamma = 1/50, len = 10, dt = 0.1, f_pdb = "pdb", 
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
        Size of the timesteps
    f_pdb: str
        Name of the pdb file (excluding the file extension, which should 
        simply be .pdb)
    f_ff: str
        Name of the forcefield file (excluding the file extension, which should 
        simply be .xml)
    """
    # set up strings
    # THE FOLLOWING FILE NAMES ARE PLACEHOLDERS: CHANGE WHEN POSSIBLE
    pdbPath = "pdb_files/{}.pdb".format(f_pdb)
    forcefieldPath = "pdb_files/{}.xml".format(f_ff)

    pdb = app.PDBFile(pdbPath)
    forcefield = app.ForceField(forcefieldPath)
    nonbonded = app.NoCutoff

    system = forcefield.createSystem(pdb.topology, nonbondedMethod = nonbonded,
                                    nonbondedCutoff = 1e3*unit.nanometer,
                                    constraints = None)
    
    beta = 1./T

    integrator = mm.LangevinIntegator(T*unit.kelvin, gamma, dt)

    platform = mm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.computeVirtualSites()
    	
    simulation.context.setVelocitiesToTemperature(Temperature*unit.kelvin)
    