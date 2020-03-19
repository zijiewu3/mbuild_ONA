import mbuild as mb
from foyer import Forcefield


class BB(mb.Compound):
    def __init__(self):
        super(BB, self).__init__(pos=[0.0, 0.0, 0.0], name='bb')
        self.r0BB = 0.084
        self.r0BBHB = 0.037
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB')
        self.add(bead)

        port = mb.Port(anchor=bead, orientation=[
                       0, 1, 0], separation=self.r0BB / 2)
        self.add(port, 'up')
        port = mb.Port(anchor=bead, orientation=[
                       0, -1, 0], separation=self.r0BB / 2)
        self.add(port, 'down')
        port = mb.Port(anchor=bead, orientation=[
                       1, 0, 0], separation=self.r0BBHB / 2)
        self.add(port, 'toHB')


class BBA(BB):
    def __init__(self):
        super(BBA, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_BBA'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_bba'


class BBT(BB):
    def __init__(self):
        super(BBT, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_BBT'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_bbt'


class BBG(BB):
    def __init__(self):
        super(BBG, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_BBG'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_bbg'


class BBC(BB):
    def __init__(self):
        super(BBC, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_BBC'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_bbc'


class HB(mb.Compound):
    def __init__(self):
        super(HB, self).__init__(pos=[0.0, 0.0, 0.0], name='hb')
        self.r0BBHB = 0.037

        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='HB')
        self.add(bead)

        port = mb.Port(
            anchor=bead, orientation=[-1, 0, 0], separation=self.r0BBHB / 2)
        self.add(port, 'toBB')


class HBA(HB):
    def __init__(self):
        super(HBA, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_HBA'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_hba'


class HBT(HB):
    def __init__(self):
        super(HBT, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_HBT'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_hbt'


class HBG(HB):
    def __init__(self):
        super(HBG, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_HBG'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_hbg'


class HBC(HB):
    def __init__(self):
        super(HBC, self).__init__()
        # Name the actual child particle that will be seen by atomtyping to '_BBA'
        for p in self.particles():
            p.name = '_HBC'
        # Name the entire compound (particle+its ports) to '_bba'
        self.name = '_hbc'


class NAT(mb.Compound):
    def __init__(self):
        super(NAT, self).__init__()
        bb = BBT()
        hb = HBT()
        self.add((bb, hb))
        mb.force_overlap(
            move_this=hb, from_positions=hb['toBB'], to_positions=bb['toHB']
        )


class NAA(mb.Compound):
    def __init__(self):
        super(NAA, self).__init__()
        bb = BBA()
        hb = HBA()

        mb.force_overlap(
            move_this=hb, from_positions=hb['toBB'], to_positions=bb['toHB']
        )

        self.add((bb, hb))


class NAG(mb.Compound):
    def __init__(self):
        super(NAG, self).__init__()
        bb = BBG()
        hb = HBG()

        self.add((bb, hb))

        mb.force_overlap(
            move_this=hb, from_positions=hb['toBB'], to_positions=bb['toHB']
        )


class NAC(mb.Compound):
    def __init__(self):
        super(NAC, self).__init__()
        bb = BBC()
        hb = HBC()

        self.add((bb, hb))

        mb.force_overlap(
            move_this=hb, from_positions=hb['toBB'], to_positions=bb['toHB']
        )


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


class DNA(mb.Compound):
    def __init__(self, sequence=None):
        super(DNA, self).__init__()
        self.sequence = sequence

        if len(self.sequence) == 0:
            raise ValueError('Sequence must not be empty')

        prev_base = get_NA(self.sequence[0])
        self.add(prev_base)

        for current_base in self.sequence[1:]:
            current_base = get_NA(current_base)
            self.add(current_base)
            # Define an convention of the down port of "earlier beads" connects
            # to the up port of "later beads" Have to move new_NA to last_NA
            # because the positions are screwed up otherwise. This is a bug and
            # should be expected to be fixed soon Have to always refer to the
            # last available port since first last_NA has two available ports
            mb.force_overlap(
                move_this=current_base,
                from_positions=(current_base.all_ports())[0],
                to_positions=(prev_base.all_ports())[-1],
            )
            prev_base = current_base

    def print_seq(self):
        return str(self.sequence)

    def visualize(self, show_ports=False):
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

        for particle in self.particles():
            if not particle.name:
                particle.name = 'UNK'

        view = py3Dmol.view()
        rad = {
            '_BBA': 0.5,
            '_BBC': 0.5,
            '_BBG': 0.5,
            '_BBT': 0.5,
            '_HBA': 0.22,
            '_HBC': 0.22,
            '_HBG': 0.22,
            '_HBT': 0.22,
        }
        col = {
            '_BBA': '0xff0000',
            '_BBC': '0x4bd1cc',
            '_BBG': '0x696969',
            '_BBT': '0xdaa520',
            '_HBA': '0x8b0000',
            '_HBC': '0x008b8b',
            '_HBG': '0x2f4f4f',
            '_HBT': '0xd2691e',
        }

        for p in self.particles(include_ports=False):
            view.addSphere(
                {
                    'center': {'x': p.pos[0], 'y': p.pos[1], 'z': p.pos[2]},
                    'radius': rad[p.name],
                    'color': col[p.name],
                    'alpha': 0.9,
                }
            )
        view.zoomTo()
        view.show()

        return view
