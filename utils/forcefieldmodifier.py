def modify_reaxff_dissociation_energy_with_ratio(filename, output_filename, energy_ratio, pbe_ratio):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Flags to identify sections in the file
    in_atom_section = False
    in_bond_section = False
    in_off_diagonal_section = False
    
    # Output lines
    modified_lines = []
    
    for line in lines:
        # Start atom section detection
        if 'Nr of atoms' in line:
            in_atom_section = True
            in_bond_section = False
            in_off_diagonal_section = True
            modified_lines.append(line)
        
            continue

        # Start bond section detection            
        if 'Nr of bonds;' in line:
            in_atom_section = False
            in_bond_section = True
            in_off_diagonal_section = False
            modified_lines.append(line)
            continue
        
        # Start off-diagonal section detection
        if 'Nr of off-diagonal terms' in line:
            in_bond_section = False  # End of bond section
            in_off_diagonal_section = True
            in_atom_section = False
            modified_lines.append(line)
            continue
        
        # End off-diagonal section detection
        if 'Nr of angles' in line:
            in_off_diagonal_section = False  # End of off-diagonal section
            in_atom_section = False
            in_bond_section = False
            modified_lines.append(line)
            continue
        
        # Process atom section (modifying Evdw values)
    
        if in_atom_section:
            # Split the line by spaces and check if it's an atom line
            pass
        # Process bond dissociation energy in bond and off-diagonal sections
        if in_bond_section or in_off_diagonal_section:
            parts = line.strip().split()
            if len(parts) == 10:  # Ensure there are enough parts
                try:
                    
                    # Keep the first two parts (integers) as they are
                    atom1 = int(parts[0])
                    atom2 = int(parts[1])
                    
                    # Modify the third part (energy) by multiplying with the energy ratio
                    current_energy1 = float(parts[2])
                    new_energy1 = current_energy1 * energy_ratio

                    current_energy2 = float(parts[3])
                    new_energy2 = current_energy2 * energy_ratio

                    current_energy3 = float(parts[4])
                    new_energy3 = current_energy3 * energy_ratio


                    pbe1 = float(parts[5])
                    new_pbe1 = pbe1 * pbe_ratio
                    
                    parts[2] = f"{new_energy1:8.4f}"
                    parts[3] = f"{new_energy2:8.4f}"
                    parts[4] = f"{new_energy3:8.4f}"
                    parts[5] = f"{new_pbe1:8.4f}"
                    formatted_line = f"{atom1:2} {atom2:2} " + " ".join(f"{float(p):8.4f}" for p in parts[2:])
                 
                    modified_lines.append(formatted_line + "\n")
                except ValueError:
                    # If parts[2] is not a valid float, leave line unchanged
                    modified_lines.append(line)
                continue

        # Add lines that are outside bond and off-diagonal sections as they are
        modified_lines.append(line)
    
    # Write the modified file
    with open(output_filename, 'w') as file:
        file.writelines(modified_lines)


import shutil

def modify_bond_parameters(file_path, factor=None, reset=False, backup_path=None):
    """
    Modify bond parameters in a reactive MD force field file.
    
    Parameters:
      file_path   : Path to the force field file.
      factor      : Multiplicative factor to apply to the bond parameters.
      reset       : If True, restore the original file from backup.
      backup_path : Path to the backup file containing default values.
    
    The function locates the bond section (after the header line 
    containing "! Nr of bonds") and updates the first three numeric 
    parameters after the atom indices.
    """
    if reset:
        if backup_path is None:
            raise ValueError("To reset, you must provide a backup_path.")
        # Overwrite file with backup contents
        shutil.copyfile(backup_path, file_path)
        return

    # Read in the original file lines.
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Locate the bonds section: the header line containing "! Nr of bonds"
    bond_header_index = None
    for i, line in enumerate(lines):
        if "! Nr of bonds" in line:
            bond_header_index = i
            break
    if bond_header_index is None:
        raise ValueError("Bond section header not found in the file.")

    # The next non-comment line should contain the number of bonds.
    num_bonds_line_index = bond_header_index + 1
    num_bonds = int(lines[num_bonds_line_index].strip().split()[0])

    # The bond parameter lines start immediately after the number-of-bonds line.
    bonds_start = num_bonds_line_index + 1
    bonds_end = bonds_start + num_bonds

    # Process each bond line.
    for j in range(bonds_start, bonds_end):
        # Split the line into tokens.
        tokens = lines[j].strip().split()
        # Assume tokens[0] and tokens[1] are atom indices;
        # tokens[2] (and tokens[3] and tokens[4] if present) are bond energies.
        for k in range(2, min(5, len(tokens))):
            try:
                # Multiply the parameter by the factor.
                new_value = float(tokens[k]) * factor
                # Format the number to match original spacing (adjust width/precision as needed).
                tokens[k] = f"{new_value:10.4f}"
            except ValueError:
                # In case conversion fails, skip updating this token.
                continue
        # Reassemble the line with a single space separator and a newline.
        lines[j] = " ".join(tokens) + "\n"

    # Write out the modified file, overwriting the original.
    with open(file_path, 'w') as f:
        f.writelines(lines)
