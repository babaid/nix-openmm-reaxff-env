from openmm import unit



def get_subset_temperature(simulation, system, atom_subset):
    # this calculates the temperature of a subset of atoms. 
    velocities = simulation.context.ggetState(getVelocities=True).getVelocities()
    subset_kinetic_energy = 0.0*unit.kilojoule/unit.mole
    for idx in atom_subset:
        mass = system.getParticleMass(idx) 
        velocity = velocities[idx]
        kinetic_energy = 0.5 * mass * (velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
        subset_kinetic_energy += kinetic_energy
    N_subset = len(subset_indices)
    k_B = 1.380649e-23 * 6.02214076e23  # J/mol/K
    temperature = (2 / (3 * k_B * N_subset)) * subset_kinetic_energy
    return temperature