#!/usr/bin/env python
# coding: utf-8

# In[7]:


from __future__ import print_function
import mbuild as mb
import builtins as __builtin__
import numpy as np
from foyer import Forcefield
#All the old import from compound
import collections
from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
import os
import sys
import tempfile
from warnings import warn
##
from mbuild.utils.io import run_from_ipython, import_
r0BB = 0.084
r0BBHB = 0.037


# In[2]:




# In[30]:




# In[3]:





# In[8]:


class BB(mb.Compound):
    def __init__(self):
        super(BB, self).__init__(pos=[0,0,0],name = 'bb')#Initizlize an instance of abstract class
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, 1, 0], separation=r0BB/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, -1, 0], separation=r0BB/2)
        self.add(port, 'down')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[1, 0, 0], separation=r0BBHB/2)
        self.add(port,'toHB')
        

class BBA(BB):
    def __init__(self):
        super(BBA,self).__init__()
        for par in self.particles():
            par.name = '_BBA'
            #print(par.name)
        
class BBT(BB):
    def __init__(self):
        super(BBT,self).__init__()
        for par in self.particles():
            par.name = '_BBT'

class BBG(BB):
    def __init__(self):
        super(BBG,self).__init__()
        for par in self.particles():
            par.name = '_BBG'
class BBC(BB):
    def __init__(self):
        super(BBC,self).__init__()
        for par in self.particles():
            par.name = '_BBC'
        
        
class HB(mb.Compound):
    def __init__(self):
        super(HB, self).__init__(pos=[0.0, 0.0, 0.0], name='hb')
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='HB')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('HB'))[0], orientation=[-1, 0, 0], separation=r0BBHB/2)
        self.add(port, 'toBB')
       

class HBA(HB):
    def __init__(self):
        super(HBA,self).__init__()
        for par in self.particles():
            par.name = '_HBA'
        
class HBT(HB):
    def __init__(self):
        super(HBT,self).__init__()
        for par in self.particles():
            par.name = '_HBT'
        
class HBG(HB):
    def __init__(self):
        super(HBG,self).__init__()
        for par in self.particles():
            par.name = '_HBG'
        
class HBC(HB):
    def __init__(self):
        super(HBC,self).__init__()
        for par in self.particles():
            par.name = '_HBC'
        
class NAT(mb.Compound):
    def __init__(self):
        super(NAT,self).__init__()
        bb = BBT()
        hb = HBT()
       # hb.translate([1,1,1])
        self.add((bb,hb))
        #print(hb['toBB'].anchor.parent)
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])

class NAA(mb.Compound):
    def __init__(self):
        super(NAA,self).__init__()
        bb = BBA()
        hb = HBA()
        #hb.translate([1,1,1])
        
        #print(hb['toBB'].anchor)
        #print(hb['toBB'].anchor.all_ports())
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])
        self.add((bb,hb))
class NAG(mb.Compound):
    def __init__(self):
        super(NAG,self).__init__()
        bb = BBG()
        hb = HBG()
        #hb.translate([1,1,1])
        self.add((bb,hb))
        #print(hb['toBB'].anchor.parent)
        #Move
        #print(hb.all_ports())
        #print(np.vstack(hb['toBB'].xyz_with_ports))
        #print(bb.all_ports())
        #print(np.vstack(bb['toHB'].xyz_with_ports))
        
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])
        
class NAC(mb.Compound):
    def __init__(self):
        super(NAC,self).__init__()
        bb = BBC()
        hb = HBC()
       # hb.translate([1,1,1])
        self.add((bb,hb))
        #print(hb['toBB'].anchor.parent)
        #Move
        #print(hb['toBB'].pos)
       # print(bb.all_ports())
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])
        

def get_NA(type='A'):
    if type == 'A':
        return NAA()
    if type == 'T':
        return NAT()
    if type == 'G':
        return NAG()
    if type == 'C':
        return NAC()
    return None

    


# In[9]:


class DNA(mb.Compound):
    sequence = None
    def __init__(self,seq=None):
        self.sequence = seq
        super(DNA,self).__init__()
        seq_len = len(seq)
        if seq_len == 0:
            pass
        #print(seq[0])
        last_NA = get_NA(seq[0])
        #last_NA['hb']
        self.add(last_NA)
        
        for letter in seq[1:]:
            new_NA = get_NA(letter)
            #new_NA.translate([1,1,1])
            #print(list(new_NA))
            self.add(new_NA)
            #Here we need to pay some attention to the sequencing
            #Define an convention of the down port of "earlier beads" connects to the up port of "later beads"
            #Have to move new_NA to last_NA because the positions are screwed up otherwise. This is a bug and should be expected to be fixed soon
            #Have to always refer to the last available port since first last_NA has two available ports
            mb.force_overlap(move_this = new_NA, from_positions=(new_NA.all_ports())[0],to_positions = (last_NA.all_ports())[-1])
            last_NA = new_NA
    def print_seq(self):
        __builtin__.print(self.sequence)
            #Questions:
            #move last NA or last_NA[0]?
            #better way to reference children?
            #Do I need to translate? How to translate?
    def ONA_visualize_py3dmol(self, show_ports=False):
        """Visualize the Compound using py3Dmol.
        Allows for visualization of a Compound within a Jupyter Notebook.
        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}
        Returns
        ------
        view : py3Dmol.view
        """
        py3Dmol = import_('py3Dmol')
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        

        for particle in self.particles():
            #particle.name = remove_digits(particle.name).upper()
            if not particle.name:
                particle.name = 'UNK'
                
        view = py3Dmol.view()
        rad = {'_BBA':.5,'_BBC':.5,'_BBG':.5,'_BBT':.5,'_HBA':0.22,'_HBC':0.22,'_HBG':0.22,'_HBT':0.22}
        col = {'_BBA':'0xff0000','_BBC':'0x4bd1cc','_BBG':'0x696969','_BBT':'0xdaa520','_HBA':'0x8b0000','_HBC':'0x008b8b'
                 ,'_HBG':'0x2f4f4f','_HBT':'0xd2691e'}
        
        for p in self.particles(include_ports=False):
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius' :rad[p.name],
                'color': col[p.name],
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view


# In[38]:


class ONA_box(mb.Compound):
    def __init__(self,sequences = [],dim = [0,0,0]):
        super(ONA_box,self).__init__()
        if len(sequences) >0 and dim[0]*dim[1]*dim[2] != len(sequences):
            dim = [len(sequences),1,1]
        seq_num = 0
        for i in range(dim[0]):
            for j in range(dim[1]):
                for k in range(dim[2]): 
                    new_DNA = DNA(sequences[seq_num])
                    new_DNA.translate_to([i*0.5,j*0.5,k*0.5])
                    self.add(new_DNA)
                    seq_num += 1
    def visualize(self):
        py3Dmol = import_('py3Dmol')
        view = py3Dmol.view()
        rad = {'_BBA':.5,'_BBC':.5,'_BBG':.5,'_BBT':.5,'_HBA':0.22,'_HBC':0.22,'_HBG':0.22,'_HBT':0.22}
        col = {'_BBA':'0xff0000','_BBC':'0x4bd1cc','_BBG':'0x696969','_BBT':'0xdaa520','_HBA':'0x8b0000','_HBC':'0x008b8b'
                 ,'_HBG':'0x2f4f4f','_HBT':'0xd2691e'}
        
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        
        for chain in self.children:
            for particle in chain.particles():
                #particle.name = remove_digits(particle.name).upper()
                if not particle.name:
                    particle.name = 'UNK'
                
        
        
                for p in chain.particles(include_ports=False):
                    view.addSphere({
                        'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                        'radius' :rad[p.name],
                        'color': col[p.name],
                        'alpha': 0.9})
        view.zoomTo()
        view.show()
        return view
    def write_lammps(self,filename):
        
        cgff = Forcefield(forcefield_files='mbuild_ONA/ONA.xml')
        #apply it to the mbuild structure
        test_box_typed =cgff.apply(self,assert_dihedral_params=False)
        #Output a LAMMPS data file
        mb.formats.lammpsdata.write_lammpsdata(test_box_typed,filename,atom_style='full')

    
    def create_lammps_input_script(self,T=4.6,script_num=1,sim_name='test.lammps',time_steps=2e6,atom_style='full',seed_1='42',seed_2='42'):
        re = import_('re')
	# Calculate dielectric constant (using equation from Ahmad's soft matter paper)
        T_real = T*(0.1/8.314)*4184
        eps_r = 249.4 - 0.788*T_real + 0.00072*(T_real**2)
        eps_r = np.round(eps_r,decimals=5)

	# Specify the name of the sample input script 
        sample_name = 'mbuild_ONA/sample.in'

	# Loop over the total number of scripts and create files 
        for i in range(script_num):
        # Specify input file name 
            filename = 'test_sample' + str(i+1) + '.in'

		# Treat the first input script separately
            if (i == 0):
			# Define substitutions
                subs = [
        				('T_ref',str(T)),
        				('eps_r',str(eps_r)),
        				('sample.dcd',sim_name + '.' + str(i+1) + '.dcd'),
        				('sample.xyz',sim_name + '.' + str(i+1) + '.xyz'),
        				('time_steps_a', str(int(time_steps + 1))),
        				('time_steps_b', str(int(time_steps))),
        				('info.dat', 'info' + str(i+1) + '.dat'),
        				('atom_style full', 'atom_style ' + atom_style),
        				('seed_1',seed_1),
        				('seed_2',seed_2)
                                        ]
            else:
# Replace relevant fields (Don't forget to replace read_data and data.lammps)
                subs = [
                    ('T_ref',str(T)),
                    ('eps_r',str(eps_r)),
                    ('read_data '+ sim_name,'read_restart restart.* remap'),
                    ('sample.dcd',sim_name + '.' + str(i+1) + '.dcd'),
                    ('sample.xyz',sim_name + '.' + str(i+1) + '.xyz'),
                    ('time_steps_a', str(int(time_steps + 1))),
                    ('time_steps_b', str(int(time_steps))),
                    ('info.dat', 'info' + str(i+1) + '.dat'),
                    ('atom_style full', 'atom_style ' + atom_style),
                    ('seed_1',seed_1),
                    ('seed_2',seed_2)
                    ]
    
		# Read sample input script 
            with open(sample_name) as fp:
                text = fp.read()
    
            lines = text.splitlines()
            for line_no,line in enumerate(lines):
                for pattern,replace in subs:
                # Do replacements
                    lines[line_no] = re.sub(pattern,replace,lines[line_no])

                # Don't re-assign velocities for i > 0 (Also, don't energy minimize)
                if (i > 0):
                    del lines[105:109]
	
                with open(filename,'w') as fp: 
                    fp.write('\n'.join(lines))
		








