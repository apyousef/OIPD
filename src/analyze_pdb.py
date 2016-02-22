import os
import sys
import copy
import optparse
import string
import numpy as np
import RTPParser as rp
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
option_parser = optparse.OptionParser()
option_parser.add_option('--id',
                        type='str',
                        help='id to attach to pdb structure')

#FIXME for now this will work, but eventually we probably want to
# make the path to the charmm27.ff directory an environment variable
# that can be set to a sensible default in a shell script.
option_parser.add_option('--charmmdir',
                        type='str',
                        help='path to charmm force field directory')

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

def _get_residue_data(charmmdir):
    aminoacids_path = os.path.join(charmmdir, 'aminoacids.rtp')
    if not os.path.exists(aminoacids_path):
        print "Error: Failed to find aminoacids.rtp file."
        sys.exit(1)
    return rp.RTP(aminoacids_path)

def _do_histidine_magic(res, rtp):
    return rtp.residues.get('HSD')

def _do_atom_magic(rtp_residue, atom):
    name = atom.name
    replacement_atom = None
    if rtp_residue.name == 'ACE':
        if name == 'H1':
            replacement_atom = rtp_residue.atoms.get('HH31')
        elif name == 'H2':
            replacement_atom = rtp_residue.atoms.get('HH32')
        elif name == 'H3':
            replacement_atom = rtp_residue.atoms.get('HH33')
    elif name == 'H':
        replacement_atom = rtp_residue.atoms.get('HN')
    elif name == 'HA3':
        replacement_atom = rtp_residue.atoms.get('HA1')
    elif name == 'HB3':
        replacement_atom = rtp_residue.atoms.get('HB1')
    elif name == 'HD3':
        replacement_atom = rtp_residue.atoms.get('HD1')
    elif name == 'HE3':
        replacement_atom = rtp_residue.atoms.get('HE1')
    elif name == 'HG3':
        replacement_atom = rtp_residue.atoms.get('HG1')
    elif rtp_residue.name == 'SER' and name == 'HG':
        replacement_atom = rtp_residue.atoms.get('HG1')
    return replacement_atom

def _calculate_center_of_charge(rtp, structure):
    total_charge = 0
    for atom in structure.get_atoms():
        name = atom.get_name()
        res = atom.get_parent()
        if res.get_resname() == 'HOH':
            continue
        rtp_residue = rtp.residues.get(res.resname)
        if not rtp_residue:
            rtp_residue = _do_histidine_magic(res, rtp)

        rtp_atom = rtp_residue.atoms.get(name)
        if not rtp_atom:
            rtp_atom = _do_atom_magic(rtp_residue, atom)

        charge = rtp_atom.charge
        print charge
        total_charge+= charge
    print "total_charge = %s" % total_charge

def translate_molecule(structure, direction):
    for atoms in structure.get_atoms():
        atoms.set_coord(atoms.get_coord() + direction)

def center_molecule(center):
    return [x * -1 for x in center]

# Bio.PDB has a matrix function that probably do the same, but it might be
# better to do it explicitly for when the linear algebra becomes more advanced, 
# and we need to optimized how things are done
def genMatrix(eAngx, eAngy, eAngz):
    TMrMatT = np.zeros((1,9))
    p = np.zeros((1,10))
    tq = np.zeros((1,4))
    a1 = 0.5 * eAngy
    a2 = 0.5 * (eAngx - eAngz)
    a3 = 0.5 * (eAngx + eAngz)
   
    q1 = np.sin(a1) * np.cos(a2)
    q2 = np.sin(a1) * np.sin(a2)
    q3 = np.cos(a1) * np.sin(a3)
    q4 = np.cos(a1) * np.cos(a3)
 
    tq[0][0] = q1
    tq[0][1] = q2
    tq[0][2] = q3
    tq[0][3] = q4
    
    k = 0
    k2 = 0
    for i in range(0,4):       
        k1 = k2
        for j in range(i,4):
            p[0][k] = 2*tq[0][k1]*tq[0][k2]
            k1 = k1 + 1
            k = k + 1
        k2 = k2 + 1

    TMrMatT[0][0] = p[0][0] + p[0][9] - 1;  TMrMatT[0][4] = p[0][4] + p[0][9] - 1;   TMrMatT[0][8] = p[0][7] + p[0][9] - 1;
    s = 1.0;    #Transpose = 1
    TMrMatT[0][1] = p[0][1] + s * p[0][8];  TMrMatT[0][3] = p[0][1] - s * p[0][8];   TMrMatT[0][2] = p[0][2] - s * p[0][6];
    TMrMatT[0][6] = p[0][2] + s * p[0][6];  TMrMatT[0][5] = p[0][5] + s * p[0][3];   TMrMatT[0][7] = p[0][5] - s * p[0][3];
    return TMrMatT.reshape((3,3))

def main():
    # validate file path
    if not args:
        print "Error: No file path provided."
        sys.exit(1)
    
    filepath1 = "/home/noel/Projects/Protein_design/Insulin/OIPD/2hiu_1H.pdb"
    filepath2 = "/home/noel/Projects/Protein_design/Insulin/OIPD/2zta_1H.pdb"
    # I added a file path for the lucine zipper molecule.
    # The program will have two structures. The first one
    # is to be centered, and the secondone is the one that
    # will be placed around the first one in many different
    # orientations. For now I will call them filepath1 and 
    # filepath2.
    filepath1 = args[0]
    if not os.path.exists(filepath1):
        print "Error: File path for molecule to be centered does not exist."
        sys.exit(1)
   
    filepath2 = args[1]
    if not os.path.exists(filepath2):
        print "Error: File path for molecule to be rotated does not exist."
        sys.exit(1)

    # validate charmm force field directory
    charmmdir = options.charmmdir
    if not charmmdir:
        print "Error: No charmm force field directory provided"
        sys.exit(1)

    if not os.path.exists(charmmdir):
        print "Error: charmm force field directory does not exist."
        sys.exit(1)

    if not os.path.isdir(charmmdir):
        print "Error: %s is not a directory" % charmmdir
        sys.exit(1)


    #TODO We should eventually perform more rigorous validation,
    # e.g. verifying file permissions.

    # assign structure id. if one is not provided, default to file name.
    structure_id = options.id
    if not structure_id:
        structure_id = os.path.splitext(os.path.basename(filepath1))[0]
    # We read both structures and place insuline at the center and 
    # leucine zipper at [5,0,0]
    structure1 = pdb_parser.get_structure('Centered One', filepath1)
    center_of_mass1 = _calculate_center_of_mass(structure1)
    #print center_of_mass1
    location_vextor1 = center_of_mass1
    translate_molecule(structure1,center_molecule(center_of_mass1))
    center_of_mass1 = _calculate_center_of_mass(structure1)
    print center_of_mass1

    structure2 = pdb_parser.get_structure('Rotated One', filepath2)
    center_of_mass2 = _calculate_center_of_mass(structure2)
    #print center_of_mass2    
    translate_molecule(structure2,center_molecule(center_of_mass2))
    location_vector2 = [45,0,0]
    translate_molecule(structure2,location_vector2)
    center_of_mass2 = _calculate_center_of_mass(structure2)
    #print center_of_mass2
    # we have two structures with possible identical chain identifiers
    # we need to make sure there are no duplicates when structure 1 and 2
    # are merged. We will relable structure2.
    # TODO: There can only be as many labels as letters in the alphabet.
    # It is unlikely that we will need more lables in the near future,
    # but constructs one day may have more chains than letters in the alphabet.
    ids = {}
    for i in string.ascii_uppercase:
        ids[i] = False
    # First We used model 0 of structure 1 and turned used_ids for that chain
    # identifier to True.
    for i in structure1[0]:
        ids[i.id] = True
    # Now we go through structure 2 and if there are any chains with the same id
    # as those found in structure 1, we will change the chain ids to something else
    for i in structure2[0]:
        if ids[i.id]:
            for j in string.ascii_uppercase:
                if not ids[j]:
                    i.id = j
                    ids[j] = True
                    break

    rtp = _get_residue_data(charmmdir)
    _calculate_center_of_charge(rtp, structure)

    # We want to join the two structures into one structure, with one model
    # and the chains of structure 1 and 2. First, deepcopy copies an object 
    # recursively
    structure3 = copy.deepcopy(structure1)
    structure3.id = 'Ensamble'
    for i in structure2.get_chains():
        structure3[0].add(i)
        
    structure_id = 0
    path_dir = '/home/noel/Projects/Protein_design/Insulin/OIPD/pdbs/'
    # DELETED: (45,45,180),(45,45,360),(45,90,180),(45,90,360),(45,135,180),(45,135,360),(45,180,180),(45,180,360)
    locations = [(45,0,0),(90,0,0),(135,0,0),(180,0,0),(225,0,0),(270,0,0),(315,0,0),(360,0,0),
                 (45,45,45),(45,45,90),(45,45,135),(45,45,225),(45,45,270),(45,45,315),
                 (45,90,45),(45,90,90),(45,90,135),(45,90,225),(45,90,270),(45,90,315),
                 (45,135,45),(45,135,90),(45,135,135),(45,135,225),(45,135,270),(45,135,315),
                 (45,180,45),(45,180,90),(45,180,135),(45,180,225),(45,180,270),(45,180,315)]
    
    for i in locations:
        RotMat = genMatrix(i[0]*np.pi/180, i[1]*np.pi/180, i[2]*np.pi/180)
        center_of_mass2 = _calculate_center_of_mass(structure2)
        translate_molecule(structure2[0],center_molecule(center_of_mass2))
        translate_molecule(structure2[0],list(np.dot(np.asarray(location_vector2),RotMat)))
        io = PDBIO()
        io.set_structure(structure3)
        io.save(path_dir+'struct_'+str(structure_id)+'.pdb')        
        #print(structure_id, center_of_mass2)
        structure_id = structure_id + 1
        
if __name__ == '__main__':
    main()
