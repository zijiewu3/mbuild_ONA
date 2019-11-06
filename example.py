#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import parmed as pmd
#from foyer import Forcefield
import mbuild as mb
from mbuild_ONA.mbuild_ONA import ONA_box
#Create topology in mbuild
import py3Dmol
import mbuild_ONA
test_box = ONA_box(['AAC','GTAC'],[1,1,2])
test_box.write_lammps('test.lammps')



# In[3]:


##Create forcefield in foyer
#cgff = Forcefield(forcefield_files='ONA.xml')
##apply it to the mbuild structure
#test_box_typed =cgff.apply(test_box,assert_dihedral_params=False)
##Output a LAMMPS data file
#mb.formats.lammpsdata.write_lammpsdata(test_box_typed,'test_output.lammps',atom_style='molecular')
