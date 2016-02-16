import os
import sys
import copy
import optparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
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

def translate_molecule(structure, direction):
    for atoms in structure.get_atoms():
        atoms.set_coord(atoms.get_coord() + direction)

def center_molecule(center):
    return [x * -1 for x in center]

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

    # We copy the structure to modify the coordinates so that
    # This second molecule is at the center of the cartesian coordinate
    # system. deepcopy copies an object recursively
    structure2 = copy.deepcopy(structure)

    translate_molecule(structure2,center_molecule(center_of_mass))
    # Confirmation that Model is centered
    center_of_mass = _calculate_center_of_mass(structure2)
    print center_of_mass 

    # This allows for the output to have two models. One for the original
    # insulin, the second one for the insulin that was centered.
    # It is practical to have two insulins or more per pdb to avoid
    # a lot of pdbs per directory when we start generating structures. It 
    # can get tricky to handle.
    structure2[0].id = 1
    structure.add(structure2[0])
    # Adding a second model to the structure has two know bugs.
    # 1. The model added with structure2[0].id = 1 is placed first on the pdb
    # 2. The model with id = 0 is place in as the second model in the pdb with id = 1
    # The output pdb does have the right manipulations, that is, it was centered, and 
    # this work for know. Let's deal with the known bugs later. 
    # I output the pdb as a test for visual inspection. 
    io = PDBIO()
    io.set_structure(structure)
    io.save('../test.pdb')
    # The next step is to load the leucine zipper in this script, and translate it
    # to a x distance from the insuline. Then generate many models of the luecine zipper
    # rotated around the insulin. I am working on that. 

if __name__ == '__main__':
    main()
