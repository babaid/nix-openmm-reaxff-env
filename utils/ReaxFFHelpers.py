import numpy as np
from openmm.app import Modeller
from openmm import *
import copy

def setup_reaxff_atoms(system: System, modeller: Modeller, force: ReaxffForce, chain_names: list[str]):
    reaxff_atoms = []
    for i, atom in zip(range(system.getNumParticles()), modeller.topology.atoms()):
        if atom.residue.chain.id in chain_names:
            force.addAtom(atom.index, atom.element.symbol, True)
            reaxff_atoms.append(atom.index)
        else:
            force.addAtom(atom.index, atom.element.symbol, False)
    return reaxff_atoms

def update_nonbonded_forces_for_reaxff(system: System, modeller: Modeller, reax_atoms: list[int], adjustment_factors=[1.0, 1.0]):
    #get nonbonded force from system
    nonbonded_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded_force = force

    #1. ALL ReaxFF-ReaxFF interactions are calculated by ReaxFF
    for i in reax_atoms:
        charge_i, sigma_i, epsilon_i = nonbonded_force.getParticleParameters(i)
        for j in reax_atoms:
            charge_j, sigma_j, epsilon_j = nonbonded_force.getParticleParameters(j)
            
            #Lorenz-Berthelot combining rule
            epsilon = np.sqrt(epsilon_i*epsilon_j)
            sigma = (sigma_i+sigma_j)/2
            
            if i!=j:
                # in theory it's still at least a hard ball
                nonbonded_force.addException(i, j, 0.0, 0.0, 0.0, True);
    
    #2. all electrostatic interactions between ReaxFF-ReaxFF/MM pairs are calculated using ReaxFF
    #3. all vdW interactions between ReaxFF-ReaxFF/MM are handled by OpenMM
    for i, atom in zip(range(system.getNumParticles()), modeller.topology.atoms()):
        charge_i, sigma_i, epsilon_i = nonbonded_force.getParticleParameters(i)
        for j in reax_atoms:
            if atom.index not in reax_atoms:
                charge, sigma_j, epsilon_j = nonbonded_force.getParticleParameters(j)
                
                #Lorenz-Berthelot combining rule
                sigma_ij = (sigma_i+sigma_j)/2
                epsilon_ij =  np.sqrt(epsilon_i*epsilon_j)
                
                adjusted_sigma =sigma_ij*adjustment_factors[0]
                adjusted_epsilon = epsilon_ij*adjustment_factors[1]
                
                #No Coulomb, possibly weakened LJ
                nonbonded_force.addException(i, j, 0.0, adjusted_sigma, adjusted_epsilon, True)
           
            
def is_in_molecule(testatoms, molatoms):
    for atom in testatoms:
        if atom in molatoms:
            return True
    return False


def remove_extra_forces(system, reax_atoms, keep_mm_forces):
    harmonic_bond_force = None
    angle_force = None
    torsion_force = None

    protein_force = HarmonicBondForce()
    protein_angle_force = HarmonicAngleForce()
    protein_torsion_force = PeriodicTorsionForce()

    # for a fully functional ReaxFF/MM simulation we need to update several forces. First off,
    # we need to remove the classical bonded MM forces from the ReaxFF atoms. 
    
    for n, force in enumerate(system.getForces()):
        if isinstance(force, HarmonicBondForce):
            harmonic_bond_force = force
            print("Number of bond forces before update: ", harmonic_bond_force.getNumBonds())
            for i in range(harmonic_bond_force.getNumBonds()):
                atom1, atom2, length, k = harmonic_bond_force.getBondParameters(i)
                if is_in_molecule([atom1, atom2], reax_atoms):
                    continue
                else:
                    protein_force.addBond(atom1, atom2, length, k)
            print("Number of bond forces after update: ", protein_force.getNumBonds())
            if keep_mm_forces:
                 harmonic_bond_force.setForceGroup(0)    
            else:
                system.removeForce(n)
            break
            
    for n, force in enumerate(system.getForces()):
        if isinstance(force, HarmonicAngleForce):
            angle_force = force
            
            print("Number of angle forces before update: ", angle_force.getNumAngles())
            for i in range(angle_force.getNumAngles()):
                atom1, atom2, atom3, angle, k = angle_force.getAngleParameters(i)
                if is_in_molecule([atom1, atom2, atom3], reax_atoms):
                    continue
                else:
                    protein_angle_force.addAngle(atom1, atom2, atom3, angle, k)
            print("Number of angle forces after update: ", protein_angle_force.getNumAngles())
            if keep_mm_forces:
                angle_force.setForceGroup(0)
            else:
                system.removeForce(n)
            break
            
    for n, force in enumerate(system.getForces()):
        if isinstance(force, PeriodicTorsionForce):
            torsion_force = force
            print("Number of torsion forces before update: ", torsion_force.getNumTorsions())
            for i in range(torsion_force.getNumTorsions()):
                atom1, atom2, atom3, atom4, periodicity, phase, k = torsion_force.getTorsionParameters(i)
                if is_in_molecule(  [atom1, atom2, atom3, atom4], reax_atoms):  
                    continue
                else:
                    protein_torsion_force.addTorsion(atom1, atom2, atom3, atom4, periodicity, phase, k)
            print("Number of torsion forces after update: ", protein_torsion_force.getNumTorsions())
            if keep_mm_forces:
                 torsion_force.setForceGroup(0)
            else:
                system.removeForce(n)
            break    
    
    protein_angle_force.setForceGroup(1)
    protein_force.setForceGroup(1)
    protein_torsion_force.setForceGroup(1)
    
    system.addForce(protein_angle_force)
    system.addForce(protein_force)
    system.addForce(protein_torsion_force)

    # following this we need to update the nonbonded force.
from openmm import NonbondedForce

def copy_nonbonded_force(original_force):
    """
    Copies an OpenMM NonbondedForce object using loops.

    Parameters:
        original_force (openmm.NonbondedForce): The original force object to copy.

    Returns:
        openmm.NonbondedForce: A new copy of the NonbondedForce.
    """
    new_force = NonbondedForce()
    
    # Copy global parameters
    new_force.setNonbondedMethod(original_force.getNonbondedMethod())
    new_force.setCutoffDistance(original_force.getCutoffDistance())
    new_force.setSwitchingDistance(original_force.getSwitchingDistance())
    new_force.setEwaldErrorTolerance(original_force.getEwaldErrorTolerance())
    new_force.setReactionFieldDielectric(original_force.getReactionFieldDielectric())
    new_force.setUseDispersionCorrection(original_force.getUseDispersionCorrection())

    # Copy per-particle parameters
    for i in range(original_force.getNumParticles()):
        charge, sigma, epsilon = original_force.getParticleParameters(i)
        new_force.addParticle(charge, sigma, epsilon)

    # Copy exception parameters
    for i in range(original_force.getNumExceptions()):
        p1, p2, chargeProd, sigma, epsilon = original_force.getExceptionParameters(i)
        new_force.addException(p1, p2, chargeProd, sigma, epsilon)

    return new_force

def get_nonbonded_forces(system, reax_atoms, atom_symbols, reaxff_force, linked_mm_atoms):

    nonbonded_force = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    reax_nonbonded_force = copy_nonbonded_force(nonbonded_force)

    for i in range(reax_nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = reax_nonbonded_force.getParticleParameters(i)
        if i in reax_atoms:
            reax_nonbonded_force.setParticleParameters(i, 0, sigma, epsilon)
            reaxff_force.addAtom(i, atom_symbols[i], 0.0, True)
        else:
            reaxff_force.addAtom(i, atom_symbols[i], charge.value_in_unit(unit.elementary_charge)*0, False)
    for i in reax_atoms:
        for j in linked_mm_atoms:
            reax_nonbonded_force.addException(i, j, 0.0, 0.0, 0.0, True)
        for j in reax_atoms:
            if i!=j:
                reax_nonbonded_force.addException(i, j, 0.0, 0.0, 0.0, True)
                
                
    reax_nonbonded_force.setForceGroup(1)
    nonbonded_force.setForceGroup(0)
    
    return reax_nonbonded_force
    
def switch_nonbonded_force(system, new_nonbonded_force):
    for i, force in enumerate(system.getForces()):
        if isinstance(force, NonbondedForce):
            system.removeForce(i)
            break
    system.addForce(new_nonbonded_force)
    
def add_allbond_reax_restraint(system, reax_atoms, atom_symbols):
    harmonic_bond_force = None
    reaxff_restraint = HarmonicBondForce()
    
    for n, force in enumerate(system.getForces()):
        if isinstance(force, HarmonicBondForce):
            harmonic_bond_force = force
            print("Number of bond forces before update: ", harmonic_bond_force.getNumBonds())
            for i in range(harmonic_bond_force.getNumBonds()):
                atom1, atom2, length, k = harmonic_bond_force.getBondParameters(i)
                if is_in_molecule([atom1, atom2], reax_atoms):
                    symbols = [atom_symbols[atom1], atom_symbols[atom2]]
                    if "H" in symbols:
                        reaxff_restraint.addBond(atom1, atom2, length, k)
                else:
                    continue
            print("Number of bond forces after update: ", reaxff_restraint.getNumBonds())
            break
    reaxff_restraint.setForceGroup(3)
    system.addForce(reaxff_restraint)