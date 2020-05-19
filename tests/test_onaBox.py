import mbuild as mb
import numpy as np
from bases_fixtures import ona_box
import pytest


def test_boxChains():
    testBox = ona_box(['ATG', 'CGTA'], dim=[2, 1, 1])
    chains = [c.sequence for c in testBox.children]
    testBox.save('onabox.gsd')
    assert set(chains) == set(['ATG', 'CGTA'])


def test_writeToLammps():
    testBox = ona_box(['ATG'], dim=[1, 1, 1])
