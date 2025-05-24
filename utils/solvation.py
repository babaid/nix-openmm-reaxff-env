import numpy as np
from openff.toolkit.topology import Molecule
from openmm import unit
from scipy.spatial.transform import Rotation
def setup_mixed_solvent(molecule: Molecule, N: int, box_size: np.ndarray, min_distance: float, randomize: bool = False):
    """
    Sets up a mixed solvent box by placing `N` copies of an OpenFF molecule in a grid or randomly within a specified box size.
    
    Parameters:
        molecule (Molecule): OpenFF molecule to place in the box.
        N (int): Number of desired copies to place.
        box_size (np.ndarray): Size of the box (x, y, z) in nanometers, as a numpy array.
        min_distance (float): Minimum distance allowed between molecules (in nanometers) to avoid clashes.
        randomize (bool): If True, molecules are placed at random positions within the box.
    
    Returns:
        List of np.ndarray: List of molecule positions, each as an (k, 3) array where k is the number of atoms.
    """
    # Create a list to hold the 3D positions of each molecule instance
    all_positions = []
    
    # Retrieve molecule's atom positions in nanometers as a numpy array
    mol_positions = molecule.to_topology().get_positions().to_openmm().value_in_unit(unit.nanometers)
    mol_positions = np.array(mol_positions)  # Convert to numpy array
    mol_size = np.ptp(mol_positions, axis=0)  # approximate molecule size along each axis

    #molecule to origin
    centroid = mol_positions.mean(axis=0)
    mol_positions-=centroid
    
    # Define grid spacing (step) based on molecule size and min_distance
    grid_step = np.max(mol_size) + min_distance
    
    # Calculate grid dimensions based on box size and grid_step
    grid_x = int(box_size[0] / grid_step)
    grid_y = int(box_size[1] / grid_step)
    grid_z = int(box_size[2] / grid_step)
    
    # Check if it's possible to fit N molecules within the grid dimensions
    max_molecules = grid_x * grid_y * grid_z
    if N > max_molecules and not randomize:
        print(f"Warning: Only {max_molecules} molecules can fit within the specified box using grid placement.")
        N = max_molecules
    
    # Function to check for clashes
    def has_clash(new_position, existing_positions, min_distance_squared):
        for pos in existing_positions:
            dist_squared = np.min(np.sum((new_position - pos) ** 2, axis=1))
            if dist_squared < min_distance_squared:
                return True
        return False
    
    # Place molecules
    min_distance_squared = min_distance ** 2  # Square distance for efficient clash checking
    count = 0
    
    if randomize:
        # Random placement within the box
        while count < N:
            # Generate a random position for the molecule's center within the box
            random_position = np.random.rand(3) * (box_size - mol_size)
            # Translate molecule's atomic positions to the random position
            new_positions = mol_positions + random_position
            
            # Check for clashes with existing molecules
            if not has_clash(new_positions, all_positions, min_distance_squared):
                # Add new molecule positions to the list if no clash is detected
                all_positions.append(new_positions)
                count += 1
    else:
        # Grid placement
        for i in range(grid_x):
            for j in range(grid_y):
                for k in range(grid_z):
                    # Calculate the grid position for the molecule's center
                    grid_position = np.array([i * grid_step, j * grid_step, k * grid_step])
                    rotation = Rotation.random().as_matrix()
                    
                    # Translate molecule's atomic positions to the grid position
                    new_positions = np.dot(mol_positions, rotation) + grid_position
                    
                    # Check for clashes with existing molecules
                    if not has_clash(new_positions, all_positions, min_distance_squared):
                        # Add new molecule positions to the list if no clash is detected
                        all_positions.append(new_positions)
                        count += 1
                        if count >= N:
                            break
                if count >= N:
                    break
            if count >= N:
                break
    if N > max_molecules:
        print(f"Warning: Only {max_molecules} molecules can fit within the specified box using grid placement.")
        N = max_molecules
    # Return final list of molecule positions
    return all_positions


import numpy as np
from scipy.spatial.transform import Rotation

import numpy as np
from scipy.spatial import KDTree

# Define the function
def place_molecules(molecules, center, protein_coords, radius, grid_spacing, min_distance):
    """
    Places molecules around a binding pocket center on a grid while avoiding clashes.

    Parameters:
    - molecules (list of np.ndarray): List of molecule coordinates as numpy arrays.
    - center (np.ndarray): 3D coordinates of the binding pocket center.
    - protein_coords (np.ndarray): Coordinates of the protein atoms.
    - radius (float): Radius around the center to place the molecules.
    - grid_spacing (float): Spacing between grid points.
    - min_distance (float): Minimum allowable distance between any two atoms.

    Returns:
    - placed_molecules (list of np.ndarray): List of placed molecule coordinates.
    """
    # Create a grid around the center
    x = np.arange(center[0] - radius, center[0] + radius, grid_spacing)
    y = np.arange(center[1] - radius, center[1] + radius, grid_spacing)
    z = np.arange(center[2] - radius, center[2] + radius, grid_spacing)
    grid_points = np.array(np.meshgrid(x, y, z)).reshape(3, -1).T

    # Initialize KDTree for protein
   
    
    protein_tree = KDTree(protein_coords)

    placed_molecules = []

    # Place the first molecule at the center
    first_molecule_coords = rotate_molecule(molecules[0]) + center
    placed_molecules.append(first_molecule_coords)

    # Build KDTree for placed molecules
   
    all_coords = protein_coords
   

    for mol_idx, mol in enumerate(molecules[1:], start=2):
        np.random.shuffle(grid_points)  # Randomize grid point order
        placed = False

        for grid_point in grid_points:
            # Apply random rotation and translate to grid point
            rotated_mol = rotate_molecule(mol) + grid_point

            # Check for clashes with existing molecules and protein
            mol_tree = KDTree(rotated_mol)
            if not mol_tree.query(all_coords, k=1)[0].min() < min_distance:
                placed_molecules.append(rotated_mol)
                all_coords = np.vstack([all_coords, rotated_mol])
                placed = True
                break

        if not placed:
            print(f"Warning: Could not place molecule {mol_idx}")

    return placed_molecules

# Helper function to apply random rotation to a molecule
def rotate_molecule(molecule):
    """
    Applies a random rotation to a molecule's coordinates.
    
    Parameters:
    - molecule (np.ndarray): Coordinates of the molecule to rotate.
    
    Returns:
    - rotated_molecule (np.ndarray): Rotated molecule coordinates.
    """
    # Generate a random rotation matrix
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.random.uniform(0, 2 * np.pi)
    psi = np.random.uniform(0, 2 * np.pi)

    Rz1 = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta),  np.cos(theta), 0],
                    [0, 0, 1]])

    Ry = np.array([[ np.cos(phi), 0, np.sin(phi)],
                   [0, 1, 0],
                   [-np.sin(phi), 0, np.cos(phi)]])

    Rz2 = np.array([[np.cos(psi), -np.sin(psi), 0],
                    [np.sin(psi),  np.cos(psi), 0],
                    [0, 0, 1]])

    rotation_matrix = Rz1 @ Ry @ Rz2

    # Apply the rotation
    rotated_molecule = molecule @ rotation_matrix.T
    return rotated_molecule



