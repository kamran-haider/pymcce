import os
import numpy as np

HEAD3_HEADER = {'FL': 0, 'occ': 1, 'crg': 2, 'Em0': 3, 'pKa0': 4, 'ne': 5,
                'nH': 6, 'vdw0': 7, 'vdw1': 8, 'tors': 9, 'epol': 10,
                'dsolv': 11, 'extra': 12, 'history': 13}


class McceEnsemble(object):
    """A class representing MCCE simulation data."""

    def __init__(self, simulation_dir, parse_configs=True):
        """
        Initializes a MCCE_ensemble object from direction where simulation data is present.

        Parameters
        ----------
        simulation_dir : string
            Path of the directory, where MCCE simulation is run.
        """
        simulation_dir = os.path.abspath(simulation_dir)
        # check if dir exists
        if not os.path.isdir(simulation_dir):
            raise IOError("%s directory does not exist", simulation_dir)

        self.simulation_dir = simulation_dir
        dir_files = os.listdir(self.simulation_dir)

        for f in dir_files:
            f = os.path.join(self.simulation_dir, f)
            #if f.endswith("head3.lst"):
            #    self.head3_data, self.conf_list = self.parse_head3lst(self, f)

    def parse_head3lst(self, head3lst):
        """
        """
        pass
        #conf_data =
        #conformer_list = []
        #with open(head3lst, "r") as h3:
        #    for line in h3.readlines()[1:]:
        #        data = line.split()
        #        conformer_list.append(data[1])
        #        conf_data = np.asarray([float(x) for x in data[2:-1]])
