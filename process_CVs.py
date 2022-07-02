import numpy as np
import pandas as pd
import os


class process_CVs:

    def __init__(self, scan_rates, mass=None, area=None, volume=None):
        self.scan_rates = scan_rates
        self.mass = mass
        self.area = area
        self.volume = volume

    def cvs_avg_and_variance(self, data):
        '''
        Calculates average potential, averager current, and current variance 
        for the second to n - 1 cycles for a given scan rate from a BioLogic 
        txt file.
        '''

        cycle_num = data['cycle number'].unique()

        current = []
        for n in cycle_num[1:-1]:
            current.append(np.asarray(
                data['<I>/mA'][data['cycle number'] == n]))
        current = pd.DataFrame(current).transpose()
        avg_current = np.mean(current, axis=1)
        current_var = np.var(current, axis=1)

        potential = []
        for n in cycle_num[1:-1]:
            potential.append(np.asarray(
                data['Ewe/V'][data['cycle number'] == n]))
        potential = pd.DataFrame(potential).transpose()
        avg_potential = potential.mean(axis=1)

        return list(avg_potential), list(avg_current), list(current_var)

    def create_and_sort_data_dict(self, file_dict):
        '''
        Creates and then sorts a dictionary of the number of cycles,
        average potentials, average currents, and current variances 
        from the files listed in the file_dict. Entries are grouped
        by scan rate.
        '''

        data_dict = {}

        for key in file_dict:
            data_dict[key] = {'cycles': [],
                              'avg_potential': [],
                              'avg_current': [],
                              'current_var': []
                              }

        for key in file_dict:

            for element in file_dict[key]:

                with open(element, 'r') as f:
                    data = pd.read_table(f)
                    cycles = max(data['cycle number']) - 2
                    avg_potential, avg_current, current_var = self.cvs_avg_and_variance(
                        data)

                data_dict[key]['cycles'].append(cycles)
                data_dict[key]['avg_potential'].append(avg_potential)
                data_dict[key]['avg_current'].append(avg_current)
                data_dict[key]['current_var'].append(current_var)

        for key in data_dict:

            for i in range(len(data_dict[key]['avg_potential'])):

                for idx, val in enumerate(data_dict[key]['avg_potential'][i]):

                    if data_dict[key]['avg_potential'][i][0] != max(data_dict[key]['avg_potential'][i]):

                        pop_potential = data_dict[key]['avg_potential'][i].pop(
                            0)
                        data_dict[key]['avg_potential'][i].append(
                            pop_potential)

                        pop_current = data_dict[key]['avg_current'][i].pop(0)
                        data_dict[key]['avg_current'][i].append(pop_current)

                        pop_current_var = data_dict[key]['current_var'][i].pop(
                            0)
                        data_dict[key]['current_var'][i].append(
                            pop_current_var)

        return data_dict

    def check_and_downsample(self, data_dict):
        '''
        Checks number of data points for each entry in the data_dict. 
        Downsamples any entry that is longer than the shortet entry. 
        '''

        lengths = []

        for key in data_dict:
            for idx, val in enumerate(data_dict[key]['avg_potential']):
                lengths.append(len(val))

        for key in data_dict:
            for idx, val in enumerate(data_dict[key]['avg_potential']):
                if len(val) != min(lengths):
                    points_to_drop = np.round(np.linspace(0, len(data_dict[key]['avg_potential'][idx]) - 1,
                                                          len(data_dict[key]['avg_potential'][idx]) -
                                                          min(lengths))).astype(int)

                    potential_update = [
                        data_dict[key]['avg_potential'][idx][i] for i in points_to_drop]
                    current_update = [data_dict[key]['avg_current'][idx][i]
                                      for i in points_to_drop]
                    var_update = [data_dict[key]['current_var'][idx][i]
                                  for i in points_to_drop]

                    potential_dropped = [item for item in data_dict[key]
                                         ['avg_potential'][idx] if item not in potential_update]
                    current_dropped = [item for item in data_dict[key]
                                       ['avg_current'][idx] if item not in current_update]
                    var_dropped = [item for item in data_dict[key]
                                   ['current_var'][idx] if item not in var_update]

                    data_dict[key]['avg_potential'][idx] = potential_dropped
                    data_dict[key]['avg_current'][idx] = current_dropped
                    data_dict[key]['current_var'][idx] = var_dropped

    def get_weighted_avgs_std(self, data_dict):
        '''
        Generates dictionary containing the weighted averages of the
        potenials and currents and the current standard deviation for
        each scan rate in the data_dict. Weighted averages use the
        number of cycles recorded at each scan rate
        '''

        final_data = {}

        for key in data_dict:

            weighted_avg_potential = np.zeros(
                len(data_dict[key]['avg_potential'][0])
            )
            weighted_avg_current = np.zeros(
                len(data_dict[key]['avg_current'][0])
            )
            averaged_current_variance = np.zeros(
                len(data_dict[key]['current_var'][0])
            )

            if self.mass != None:
                weighted_avg_current_density = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                averaged_current_density_variance = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_spec_capacitance = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_spec_capacitance_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_spec_capacity = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_spec_capacity_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_spec_charge = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_spec_charge_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )

            if self.area != None:
                weighted_avg_areal_capacitance = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_areal_capacitance_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_areal_capacity = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_areal_capacity_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_areal_charge = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_areal_charge_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )

            if self.volume != None:
                weighted_avg_volume_capacitance = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_volume_capacitance_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_volume_capacity = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_volume_capacity_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )
                weighted_avg_volume_charge = np.zeros(
                    len(data_dict[key]['avg_current'][0])
                )
                avg_volume_charge_std_dev = np.zeros(
                    len(data_dict[key]['current_var'][0])
                )

            for idx, val in enumerate(data_dict[key]['avg_potential']):
                weighted_avg_potential += (np.asarray(val) * (
                    data_dict[key]['cycles'][idx] / sum(data_dict[key]['cycles'])))

            for idx, val in enumerate(data_dict[key]['avg_current']):
                weighted_avg_current += (np.asarray(val) * (
                    data_dict[key]['cycles'][idx] / sum(data_dict[key]['cycles'])))

            for idx, val in enumerate(data_dict[key]['current_var']):
                averaged_current_variance += np.asarray(val)

            final_data[key] = {
                'w_avg_potential': weighted_avg_potential,
                'w_avg_current': weighted_avg_current,
                'current_std_dev': np.sqrt(averaged_current_variance)
            }

        return final_data
