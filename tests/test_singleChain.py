import mbuild as mb
import numpy as np
import pytest


def test_sc_sequence():
    chain = mb.recipes.DNA('ATGC')
    name_list = []
    for na in chain.children:
        name_list.append(na.name)
    assert name_list == ['NAA', 'NAT', 'NAG',
                         'NAC'] and chain.print_seq() == 'ATGC'


def test_sc_particles():
    chain = mb.recipes.DNA('ATGC')
    name_list = []
    for na in chain.children:
        for p in na.children:
            name_list.append(p.name)
    assert name_list == ['_bba', '_hba', '_bbt',
                         '_hbt', '_bbg', '_hbg', '_bbc', '_hbc']


def test_HBBB_bond():
    # Check there is only 1 bond that has the only two particles in the system connected
    chain = mb.recipes.DNA('A')
    particle_set = set([p for p in chain.particles()])
    bonds = [b for b in chain.bonds()]
    firstBond = set(bonds[0])
    assert len(bonds) == 1 and firstBond == particle_set


def test_bonds_twoResidue():
    chain = mb.recipes.DNA('AG')
    pl = [p for p in chain.particles()]
    bonds = [b for b in chain.bonds()]
    firstBond = frozenset(bonds[0])
    secondBond = frozenset(bonds[1])
    thirdBond = frozenset(bonds[2])
    assert len(bonds) == 3 and set([firstBond, secondBond, thirdBond]) == set(
        [frozenset([pl[0], pl[1]]), frozenset([pl[2], pl[3]]), frozenset([pl[0], pl[2]])])
