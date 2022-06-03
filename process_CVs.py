import numpy as np
import pandas as pd
import os


class process_CVs:

    def __init__(self, rootdir, scan_rates):
        self.rootdir = rootdir
        self.scan_rates = scan_rates

    def create_file_dict(self):
        ''''''

        file_dict = {}

        directory_path = [x[0] for x in os.walk(self.rootdir)][1:]
        file_path = [x[2] for x in os.walk(self.rootdir)][1:]

        full_path_list = []
        for idx, nested_list in enumerate(file_path):
            temp = []
            for element in nested_list:
                temp.append('{}\\{}'.format(directory_path[idx], element))
            full_path_list.append(temp)

        for scan_rate, path in zip(self.scan_rates, full_path_list):
            file_dict[scan_rate] = path

        return file_dict

    def cvs_avg_and_variance(self, data):
        '''Calculates average potential, current, and current variance from all cycles
        (excluding first and last) for a given scan rate from a BioLogic txt file.'''

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
        '''    '''

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
        ''''''

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

        final_data = {}

        for key in data_dict:

            weighted_avg_potential = np.zeros(
                len(data_dict[key]['avg_potential'][0])
            )
            weighted_avg_current = np.zeros(
                len(data_dict[key]['avg_current'][0])
            )
            avg_current_std = np.zeros(
                len(data_dict[key]['current_var'][0])
            )

            for idx, val in enumerate(data_dict[key]['avg_potential']):
                weighted_avg_potential += (
                    np.asarray(val) * (
                        data_dict[key]['cycles'][idx] /
                        sum(data_dict[key]['cycles'])
                    )
                )

            for idx, val in enumerate(data_dict[key]['avg_current']):
                weighted_avg_current += (
                    np.asarray(val) * (
                        data_dict[key]['cycles'][idx] /
                        sum(data_dict[key]['cycles'])
                    )
                )

            for idx, val in enumerate(data_dict[key]['current_var']):
                avg_current_std += np.asarray(val)

            final_data[key] = {
                'w_avg_potential': weighted_avg_potential,
                'w_avg_current': weighted_avg_current,
                'current_std_dev': np.sqrt(avg_current_std)
            }

        return final_data
