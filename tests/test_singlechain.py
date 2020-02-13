import mbuild as mb
import numpy as np
import pytest

def test_sc_sequence():
    chain = mb.recipes.DNA('ATGC')
    name_list = []
    for na in chain.children:
        name_list.append(na.name)
    assert name_list == ['NAA','NAT','NAG','NAC']

def test_sc_particles():
    chain = mb.recipes.DNA('ATGC')
    name_list = []
    for na in chain.children:
        for p in na.children:
            name_list.append(p.name)
    print(name_list[0])
    print(name_list[1])
    assert name_list == ['_BBA','_HBA','_BBT','_HBT','_BBG','_HBG','_BBC','_HBC']

def test_HBBB_bond():
    chain = mb.recipes.DNA('A')
    particle_list = []
    for p in chain.particles():
        particle_list.append(p)
    bond_gen = chain.bonds()
    firstbond = next(bond_gen)
    assert firstbond[1] == particle_list[0] and firstbond[0] == particle_list[1]
