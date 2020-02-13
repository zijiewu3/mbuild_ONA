import mbuild as mb
import numpy as np
from foyer import Forcefield

from ONA_box import ONA_box
class ONA_box(mb.Compound):
    def __init__(self,sequences = [],dim = [0,0,0]):
        super(ONA_box,self).__init__()
        if len(sequences) >0 and dim[0]*dim[1]*dim[2] != len(sequences):
            dim = [len(sequences),1,1]
        seq_num = 0
        for i in range(dim[0]):
            for j in range(dim[1]):
                for k in range(dim[2]): 
                    new_DNA = mb.recipes.DNA(sequences[seq_num])
                    new_DNA.translate_to([i*0.5,j*0.5,k*0.5])
                    self.add(new_DNA)
                    seq_num += 1
testBox = ONA_box(['ATAT','CGCG'],[2,1,1]) 
testBox.visualize()
#cgff = Forcefield(forcefield_files='ONA.xml')
        #apply it to the mbuild structure
#test_box_typed =cgff.apply(testBox,assert_dihedral_params=False)
        #Output a LAMMPS data file
#mb.formats.lammpsdata.write_lammpsdata(test_box_typed,'test.lammpsdata',atom_style='full')


