import os
import sys
import optparse
from Bio.PDB.PDBParser import PDBParser

option_parser = optparse.OptionParser()
option_parser.add_option('--id',
                        type='str',
                        help='id to attach to pdb structure')

options, args = option_parser.parse_args()

#TODO PDBParser has additional keyword arguments like PERMISSIVE.
# We may decide to add options to the option parser for such
# keywords, but for now let's not over-engineer.
pdb_parser = PDBParser()

def _calculate_center_of_mass(structure):
    total_mass = 0
    mx_total = 0
    my_total = 0
    mz_total = 0
    for atom in structure.get_atoms():
        coords = atom.coord.tolist()
        mass = atom.mass
        total_mass += mass
        mx_total += coords[0] * mass
        my_total += coords[1] * mass
        mz_total += coords[2] * mass
    return [mx_total/total_mass, my_total/total_mass, mz_total/total_mass]

def main():
    # validate file path
    if not args:
        print "Error: No file path provided."
        sys.exit(1)

    filepath = args[0]
    if not os.path.exists(filepath):
        print "Error: File path does not exist."
        sys.exit(1)

    #TODO We should eventually perform more rigorous validation,
    # e.g. verifying file permissions.

    # assign structure id. if one is not provided, default to file name.
    structure_id = options.id
    if not structure_id:
        structure_id = os.path.splitext(os.path.basename(filepath))[0]

    structure = pdb_parser.get_structure(structure_id, filepath)
    center_of_mass = _calculate_center_of_mass(structure)
    print center_of_mass

if __name__ == '__main__':
    main()