import mbuild as mb
import numpy as np


def test_box():
    testBox = mb.Box(mins=[0, 0, 0], maxs=[4, 4, 4])
    assert testBox.lengths.all() == np.array([4, 4, 4]).all()
