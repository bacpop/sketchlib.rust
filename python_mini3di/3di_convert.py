
import mini3di
from Bio.PDB import PDBParser

encoder = mini3di.Encoder()
parser = PDBParser(QUIET=True)

def pdb_to_3di(struct_name: str, filename: str) -> str:
    struct = parser.get_structure(struct_name, filename)

    structure_string = ""
    for chain in struct.get_chains():
        states = encoder.encode_chain(chain)
        sequence = encoder.build_sequence(states)
        structure_string += f"{sequence},"
    structure_string.removesuffix(",")

    return structure_string