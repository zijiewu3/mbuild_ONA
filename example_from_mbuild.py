import mbuild as mb
test_box = mb.recipes.ONA_box(['AATT','GGCC'],[2,1,1])
test_box.write_lammps('test_from_mb.lammps')
