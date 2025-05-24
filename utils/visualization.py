from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from IPython.display import SVG

def visualize_rdkit(rdkit_molecule, molecule, molSize=(300, 300), kekulize=True):
    # Create a copy of the molecule to avoid altering the original
    mc = Chem.Mol(rdkit_molecule.ToBinary())
    
    # Try to kekulize the molecule for better visual representation
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except Exception as e:
            print(f"Kekulization failed: {e}")
            mc = Chem.Mol(rdkit_molecule.ToBinary())
    
    # Ensure the molecule has 2D coordinates
    rdDepictor.Compute2DCoords(mc)
    
    # Map custom atom labels from `molecule`
    for atom1, atom2 in zip(molecule.atoms, mc.GetAtoms()):
        atom2.SetProp('atomLabel', atom1.name)
    
    # Draw the molecule in SVG format
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    
    # Clean up the SVG text for rendering
    svg = drawer.GetDrawingText()
    return SVG(svg.replace('svg:', ''))