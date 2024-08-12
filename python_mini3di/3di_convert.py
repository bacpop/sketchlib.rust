import os, sys
import mini3di
from Bio.PDB import PDBParser
from warnings import warn


encoder = mini3di.Encoder()
parser = PDBParser(QUIET=True)
#parser = PDBParser()

def pdb_to_3di(struct_name: str, filename: str) -> str:
    struct = parser.get_structure(struct_name, filename)

    #print("filename: " + filename + " structname: " + struct_name)

    structure_string = ""
    for chain in struct.get_chains():
        #if len(list(chain.get_residues())) < 3: continue
        #print(len(list(chain.get_residues())))
        try:
            states = encoder.encode_chain(chain)
            sequence = encoder.build_sequence(states)
            structure_string += f"{sequence},"
        except IndexError:
            warn("Not able to code into 3Di chain {} from protein ID {}".format(chain.__repr__(), struct_name), RuntimeWarning)
            continue
    structure_string = structure_string.removesuffix(",")
    #print(structure_string)
    return structure_string


def main():
    thel = sys.argv
    if len(thel) < 2:
        raise RuntimeError("FATAL: no input list of files provided!")

    if not os.path.isfile(thel[1]):
        raise RuntimeError("FATAL: not provided a valid file!!!")

    with open(thel[1], "r") as inf:
        for line in inf:
            tmpl = line.strip().split("\t")
            #print(tmpl[1])
            #print(pdb_to_3di(tmpl[0], tmpl[1]))
            pdb_to_3di(tmpl[0], tmpl[1])


if __name__ == '__main__':
    main()

