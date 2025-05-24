import numpy as np

def check_clash(molecule_coords, all_coords, min_distance):
    """Return True if any atom in molecule_coords is closer than min_distance to any atom in all_coords."""
    for atom in molecule_coords:
        for other in all_coords:
            if np.linalg.norm(atom - other) < min_distance:
                return True
    return False

def is_cylinder_clear(new_mol_coords, base_center, protein_coords, clearance):
    """
    Check if a cylindrical region between the new molecule's center and the base fragment center is free of protein.
    
    The cylinder has:
      - Its axis from the new molecule's center (computed as the mean of new_mol_coords) to base_center.
      - A radius given by the clearance parameter.
    
    Parameters:
    - new_mol_coords (np.ndarray): Coordinates of the new fragment.
    - base_center (np.ndarray): Center of the base fragment.
    - protein_coords (np.ndarray): Protein atom coordinates.
    - clearance (float): Minimum allowed distance from the cylinder axis (i.e. cylinder radius).
    
    Returns:
    - bool: True if the entire cylinder is free of protein atoms, False otherwise.
    """
    # Calculate the center of the new molecule
    new_center = np.mean(new_mol_coords, axis=0)
    axis_vector = base_center - new_center
    axis_length = np.linalg.norm(axis_vector)
    if axis_length == 0:
        return True
    # Normalize the axis
    axis_unit = axis_vector / axis_length

    # Loop over each protein atom and check if it falls inside the cylinder.
    for prot_atom in protein_coords:
        # Vector from the new molecule center to the protein atom
        w = prot_atom - new_center
        # Projection of w onto the axis
        t = np.dot(w, axis_unit)
        # Only consider protein atoms between the two centers
        if 0 <= t <= axis_length:
            # Compute the perpendicular distance from the protein atom to the axis
            projection = new_center + t * axis_unit
            perpendicular_distance = np.linalg.norm(prot_atom - projection)
            if perpendicular_distance < clearance:
                return False
    return True

def place_molecules(molecules, center, protein_coords, radius, grid_spacing, min_distance):
    """
    Places molecules around a binding pocket center on a grid while avoiding clashes and ensuring that
    each placed molecule has a clear cylindrical path (free of protein atoms) to the base fragment.

    Parameters:
    - molecules (list of np.ndarray): List of molecule coordinates.
    - center (np.ndarray): 3D coordinates of the binding pocket center (base fragment position).
    - protein_coords (np.ndarray): Coordinates of the protein atoms.
    - radius (float): Radius around the center to place the molecules.
    - grid_spacing (float): Spacing between grid points.
    - min_distance (float): Minimum allowable distance between any two atoms (also used as clearance for the cylinder).

    Returns:
    - placed_molecules (list of np.ndarray): List of placed molecule coordinates.
    """
    # Create a grid around the center
    x = np.arange(center[0] - radius, center[0] + radius, grid_spacing)
    y = np.arange(center[1] - radius, center[1] + radius, grid_spacing)
    z = np.arange(center[2] - radius, center[2] + radius, grid_spacing)
    grid_points = np.array(np.meshgrid(x, y, z)).reshape(3, -1).T

    placed_molecules = []
    
    # Place the first molecule (base fragment) at the center
    first_molecule_coords = rotate_molecule(molecules[0]) + center
    placed_molecules.append(first_molecule_coords)
    
    # All occupied space: protein atoms plus the base fragment's atoms
    all_coords = protein_coords.copy()
    all_coords = np.vstack([all_coords, first_molecule_coords])

    # For each remaining molecule, try to place it at a grid point.
    for mol_idx, mol in enumerate(molecules[1:], start=2):
        np.random.shuffle(grid_points)  # Randomize grid point order
        placed = False

        for grid_point in grid_points:
            # Apply random rotation and translate to grid point
            rotated_mol = rotate_molecule(mol) + grid_point

            # Check for clashes with already placed atoms
            if check_clash(rotated_mol, all_coords, min_distance):
                continue  # Clash found; skip this grid point

            # Check that the cylinder from the new molecule's center to the base center is free of protein atoms
            if not is_cylinder_clear(rotated_mol, center, protein_coords, min_distance):
                continue  # Protein is obstructing the cylindrical path

            # If both checks pass, accept the placement.
            placed_molecules.append(rotated_mol)
            all_coords = np.vstack([all_coords, rotated_mol])
            placed = True
            break

        if not placed:
            print(f"Warning: Could not place molecule {mol_idx}")

    return placed_molecules

# Dummy implementation of rotate_molecule for demonstration.
def rotate_molecule(mol_coords):
    # In practice, you might apply a random rotation here.
    # For this example, we'll just return the original coordinates.
    return mol_coords
