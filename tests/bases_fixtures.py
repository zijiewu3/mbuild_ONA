#!/usr/bin/env python
# coding: utf-8

# In[7]:


from __future__ import print_function
import mbuild as mb
import builtins as __builtin__
import numpy as np
from foyer import Forcefield
# All the old import from compound
# import collections
# from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
import pytest
import os
import sys
import tempfile
from warnings import warn
##
from mbuild.utils.io import run_from_ipython, import_
import pytest


def ona_box(sequence, dim):
    class ONA_box(mb.Compound):
        def __init__(self, sequences=[], dim=[0, 0, 0]):
            super(ONA_box, self).__init__()
            if len(sequences) > 0 and dim[0] * dim[1] * dim[2] != len(sequences):
                dim = [len(sequences), 1, 1]
                print(
                    'Number of chains and dimension do not match. Reconfiguring dim to ' + str(dim))
            seq_num = 0
            for i in range(dim[0]):
                for j in range(dim[1]):
                    for k in range(dim[2]):
                        new_DNA = mb.recipes.DNA(sequences[seq_num])
                        new_DNA.translate_to([i * 0.5, j * 0.5, k * 0.5])
                        self.add(new_DNA)
                        seq_num += 1

        def visualize(self):
            py3Dmol = import_('py3Dmol')
            view = py3Dmol.view()
            rad = {'_BBA': 0.05, '_BBC': 0.05, '_BBG': 0.05, '_BBT': 0.05,
                   '_HBA': 0.022, '_HBC': 0.022, '_HBG': 0.022, '_HBT': 0.022}
            col = {'_BBA': '0xff0000', '_BBC': '0x4bd1cc', '_BBG': '0x696969', '_BBT': '0xdaa520',
                   '_HBA': '0x8b0000', '_HBC': '0x008b8b', '_HBG': '0x2f4f4f', '_HBT': '0xd2691e'}

            def remove_digits(x): return ''.join(
                i for i in x if not i.isdigit() or i == '_')

            # DOwn to chain
            for chain in self.children:
                # Down to nucleotic acid
                for NA in chain.children:
                    # Down to particles
                    for p in NA.children:
                        print(p.name)
                        view.addSphere({
                            'center': {'x': p.pos[0], 'y': p.pos[1], 'z': p.pos[2]},
                            'radius': rad[p.name],
                            'color': col[p.name],
                            'alpha': 0.9})
            view.zoomTo()
            view.show()
            return view

        def write_lammps(self, filename):

            cgff = Forcefield(forcefield_files='ONA.xml')
            # apply it to the mbuild structure
            test_box_typed = cgff.apply(self, assert_dihedral_params=False)
            # Output a LAMMPS data file
            mb.formats.lammpsdata.write_lammpsdata(
                test_box_typed, filename, atom_style='full')
    return ONA_box(sequence, dim)
