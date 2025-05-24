#These are tools to handle some things for topologies, like naming for example
from openmm.app import Topology
from openmm import Vec3
class TopologyTools:
    @staticmethod
    def add_chain_name(chain_id:str, topology: Topology):
        for i, chain in enumerate(topology.chains()):
            if chain.id != "X":
                print("Warning, renaming chain, old chain id was ", chain.id, ".")
            if i>0:
                chain.id = f"{chain_id}-{i}"
            else:
                chain.id = chain_id
    @staticmethod
    def get_atom_index(topology: Topology, chain_id, atom_name: str):
        for chain in topology.chains():
            if chain.id == chain_id:
                for atom in chain.atoms():
                    if atom.name == atom_name:
                        return atom.index


def calculate_bounding_box(pdb):
    positions = pdb.positions
    
    min_x = min(pos.x for pos in positions)
    min_y = min(pos.y for pos in positions)
    min_z = min(pos.z for pos in positions)
    
    max_x = max(pos.x for pos in positions)
    max_y = max(pos.y for pos in positions)
    max_z = max(pos.z for pos in positions)
    
    return Vec3(min_x, min_y, min_z), Vec3(max_x, max_y, max_z)
def save_system_state(system, simulation, modeller, output, identifier):

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(os.path.join(output, 'systems', f'system_{identifier}.xml'), 'w') as outputfile:
        outputfile.write(XmlSerializer.serialize(state))
    
    outputfile = os.path.join(output, 'topologies', f'topology_{identifier}.pdb')
    PDBFile.writeFile(modeller.topology, state.getPositions(), outputfile)
    
    with open(os.path.join(output, 'states', f'state_{identifier}.xml'), 'w') as outputfile:
        outputfile.write(XmlSerializer.serialize(system))
