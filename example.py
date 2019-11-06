import mbuild as mb
from mbuild_ONA.mbuild_ONA import ONA_box
import py3Dmol
import mbuild_ONA
test_box = ONA_box(['AAC','GTAC'],[1,1,2])
test_box.write_lammps('test.lammps')


