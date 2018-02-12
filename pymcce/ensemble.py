import os
import numpy as np

HEAD3_HEADER = {'FL': 0, 'occ': 1, 'crg': 2, 'Em0': 3, 'pKa0': 4, 'ne': 5,
                'nH': 6, 'vdw0': 7, 'vdw1': 8, 'tors': 9, 'epol': 10,
                'dsolv': 11, 'extra': 12, 'history': 13}


class McceEnsemble(object):
    """A class representing MCCE simulation data."""

    def __init__(self, ms_data_file, head3lst_file, fort38_file):
        self.ms_data_file = ms_data_file
        self.byte_indices = None
        self.total_microstates = 0
        self.total_records = 0
        res_list = []

        with open(self.ms_data_file, "rb") as md:
            bytes_n_res = md.read(4)
            n_res = struct.unpack('i', bytes_n_res)[0]
            for i in range(n_res):
                resname = str(md.read(8))
                res_list.append(resname)
            self.n_res = n_res
        self.residue_list = res_list
        # self.residue_hb_matrix = np.zeros((n_res, n_res), dtype=float)

        # with open(ms_gold_file, "r") as ms_gold:
        #    self.n_res = len([res.strip() for res in ms_gold.readlines() if len(res) != 0])
        conf_data = {}
        conf_id_name_map = {}
        with open(head3lst_file, "r") as h3:
            for line in h3.readlines()[1:]:
                data = line.split()
                conf_id = data[1]
                conf_data[conf_id] = [int(data[0]), 0.0, []]
                conf_id_name_map[int(data[0])] = conf_id
        self.conformer_data = conf_data
        self.conf_id_name_map = conf_id_name_map

        # read fort38 data
        with open(fort38_file, 'r') as f:
            lines = [l.strip().split() for l in f.readlines()]
            header = lines.pop(0)
            for index, l in enumerate(lines):
                conf_id = index + 1
                conf_name = l[0]
                occ = float(l[1])
                if conf_name in self.conformer_data.keys():
                    self.conformer_data[conf_name][1] += occ

        residue_data = {}
        for k in self.conformer_data.keys():
            conf_name = k
            res_key = conf_name[:3] + conf_name[5:10]
            if "CTR" not in res_key and "NTR" not in res_key:
                if res_key not in residue_data.keys():
                    residue_data[res_key] = [(k, self.conformer_data[k][0])]
                else:
                    residue_data[res_key].append((k, self.conformer_data[k][0]))

        self.residue_data = residue_data

    def generate_byte_indices(self, sample_frequency=100, filename=None):
        """
        Generate byte indices

        Parameters
        ----------
        n_res : TYPE
            Description
        sample_frequency : int, optional
            Description
        filename : None, optional
            Description

        Returns
        -------
        rec_indices : list
            A list of the starting bytes for each record in the microstate data file
        """
        if filename is None:
            filename = self.ms_data_file

        start_byte = 4 + (8 * self.n_res)
        bytes_per_record = (self.n_res * 2) + 20
        file_size = os.path.getsize(filename)
        # n_records = (file_size - start_byte) / bytes_per_record
        rec_indices = list(range(start_byte, file_size, sample_frequency * bytes_per_record))
        self.total_records = len(rec_indices)
        self.byte_indices = rec_indices

    def parse_records(self):
        """
        Parse ms.dat
        """
        trajectory = np.zeros([self.total_records, self.n_res], dtype=int)
        state_counts = np.zeros([self.total_records], dtype=int)
        energies = np.zeros([self.total_records], dtype=float)
        # print trajectory
        progress_counter = 0
        # print_progress_bar(progress_counter, self.total_records)
        with open(self.ms_data_file, "rb") as ms:
            for index, record in enumerate(self.byte_indices):
                ms.seek(record)
                bytes_conf_ids = ms.read(2 * self.n_res)
                bytes_energies_1 = ms.read(8)
                ms.seek(ms.tell() + 8)
                energy = struct.unpack("d", bytes_energies_1)[0]
                bytes_state_count = ms.read(4)
                trajectory[index, :] = np.asarray(struct.unpack(str(self.n_res) + "H", bytes_conf_ids))
                # print(struct.unpack(str(self.n_res) + "H", bytes_conf_ids)[-2:])
                state_count = struct.unpack("i", bytes_state_count)[0]
                self.total_microstates += state_count
                state_counts[index] += state_count
                energies[index] += energy
                progress_counter += 1
                # print_progress_bar(progress_counter, self.total_records)
        self.trajectory = trajectory
        self.state_counts = state_counts
        self.energies = energies

    def parse_struct(self, step2out_pdb, n_conf):
        """
        Obtain structural data for each conformer
        """
        print("Parsing structure ...")
        st = md.load_pdb(step2out_pdb, frame=0, no_boxchk=True)
        for index, res in enumerate(st.topology.residues):
            if res.name == "HOH":
                k = "%sA%04d" % (res.name, res.resSeq)
                residue_at_ids = [at.index for at in res.atoms]
                for c in range(n_conf):
                    conf_name = "HOH01W%04d_%03d" % (res.resSeq, c + 1)
                    conf_rows = c * 3, (c * 3) + 3
                    conf_at_ids = residue_at_ids[conf_rows[0]:conf_rows[1]]
                    # conf_st = st.atom_slice(conf_at_ids)
                    self.conformer_data[conf_name][2].extend(conf_at_ids)
                    # angle = md.compute_angles(st, [[conf_at_ids[1], conf_at_ids[0], conf_at_ids[2]]])
                    # angle = md.utils.in_units_of(angle, "radians", "degrees")[0]
                    # if angle <= 104:
                    # print(conf_name, angle, conf_at_ids)
        self.structure = st
        print("Done.")

    def calculate_water_dipoles(self):
        charge_dict = {"O": -2.000, "1H": 1.000, "2H": 1.000}
        chg = [-0.834, 0.417, 0.417]
        # chg = [-2.0, 1.0, 1.0]

        v = np.array([0, 0, 1])
        for index, c in enumerate(self.conformer_data.keys()):
            if len(self.conformer_data[c][-1]) != 0:
                conf_ids = self.conformer_data[c][-1]
                pos = self.structure.xyz[0, conf_ids, :].T
                pos *= -10.0
                dp = pos.dot(chg)
                # calculate projection
                proj = np.multiply(dp.dot(v) / v.dot(v), v)
                scaling_factor = dp.dot(v) / v.dot(v)
                # print(index, scaling_factor)
                # snapshot_st = msa.structure.atom_slice(conf_ids)
                # snapshot_st.save_pdb("%d_wat_conf.pdb" % index)

                self.conformer_data[c].append(scaling_factor)
            else:
                self.conformer_data[c].append(None)

    def parse_energy_data(self):

        energies_dir = "/Users/kamranhaider/Dropbox/ProtonHopping_scratch/data/mcce_water_profiling/5_conf_per_water/energies"
        opp_files = [f for f in os.listdir(energies_dir) if f.endswith(".opp")]
        opp_files = sorted(opp_files)
        n = len(self.conformer_data.keys())
        print(len(self.conformer_data.keys()), len(opp_files))
        pairwise_energies = np.zeros((n, n))
        for c in self.conformer_data.keys():
            # opp = os.path.join(energies_dir, f)
            print(c, opp_files[0], opp_files[1])
            break
            # with open(opp, "r") as o:
            #    data = o.readlines()
