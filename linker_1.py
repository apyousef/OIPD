# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 00:16:21 2016

@author: noel
"""

from Bio.PDB import *

parser = PDBParser()

structure = parser.get_structure('Insulin', '/home/noel/Projects/Protein_design/Insulin/OIPD/2hiu_1H.pdb')

''' TO DO:
    1. Create a class to read charmm27.ff/ffnonbonded.itp for masses. Place in a dictionary of atom types.
    2. Pass dicctionary with atom types and structure[0] to return center of mass from structure_values()
    3. Create a new structure so that both the original structure plus the centered (by center of mass)
       is output in the same pdb.
'''

import centers as ct
newCenter = structure_values()
CM = [1,1,1] #newCenter.center_of_mass(structure[0])

print("Done calling centers.py", CM)

for atom in structure.get_atoms():
    print(atom.get_coord(), CM) 
    atom.set_coord(atom.get_coord() - CM)
    print(atom.get_coord())

#eAngx = 0.0;      eAngy = 0.0;       eAngz = 0.0;
#a1 = 0.5 * eAngy;
#a2 = 0.5 * (eAngx - eAngz);
#a3 = 0.5 * (eAngx + eAngz);
#System.out.println("   a123 = "+a1+" "+a2+" "+a3);
#TM.q_u1[0] = Math.sin(a1) * Math.cos(a2);
#TM.q_u2[0] = Math.sin(a1) * Math.sin(a2);
#TM.q_u3[0] = Math.cos(a1) * Math.sin(a3);
#TM.q_u4[0] = Math.cos(a1) * Math.cos(a3);
#System.out.println("   q1234 = "+TM.q_u1[0]+" "+TM.q_u2[0]+" "+TM.q_u3[0]+" "+TM.q_u4[0]);
#tq[0] = TM.q_u1[0];  tq[1] = TM.q_u2[0];  tq[2] = TM.q_u3[0];    tq[3] = TM.q_u4[0];
#for(k = 0, k2 = 0; k2 < 4; k2++){
#    for(k1 = k2; k1 < 4; k1++, k++){
#        p[k] = 2.0*tq[k1]*tq[k2];
#    }
#}
#TM.rMatT[0] = p[0] + p[9] - 1;  TM.rMatT[4] = p[4] + p[9] - 1;   TM.rMatT[8] = p[7] + p[9] - 1;
#s = 1.0;    //Transpose = 1
#TM.rMatT[1] = p[1] + s * p[8];   TM.rMatT[3] = p[1] - s * p[8];   TM.rMatT[2] = p[2] - s * p[6];
#TM.rMatT[6] = p[2] + s * p[6];   TM.rMatT[5] = p[5] + s * p[3];   TM.rMatT[7] = p[5] - s * p[3];
