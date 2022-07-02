import numpy as np
import pandas as pd
import os


class rate_performance:

    def __init__(self, scan_rates, mass=None, area=None, volume=None):
        self.scan_rates = scan_rates
        self.mass = mass
        self.area = area
        self.volume = volume

    def rate_performance(self, data):
        '''
        Calculates average capacitance, capacity, and charge (and error values) 
        for a series of scan rates.Returns a plot of average capacitance with 
        error bars versus scan rate.
        '''
        potential_window = 0

        if potential_window == 0:
            potential_window += (max(data['Ewe/V']) - min(data['Ewe/V']))

        cycle_number = data['cycle number'].unique()

        charge_list = []
        charge_data = data[data['<I>/mA'] > 0]
        for n in cycle_number[1:-1]:
            charge_list.append(
                np.trapz(
                    charge_data['<I>/mA'][charge_data['cycle number'] == n],
                    charge_data['time/s'][charge_data['cycle number'] == n]
                )
            )

        charge = np.mean(charge_list)
        charge_var = np.var(charge_list)

        discharge_list = []
        discharge_data = data[data['<I>/mA'] < 0]
        for n in cycle_number[1:-1]:
            discharge_list.append(
                np.trapz(
                    discharge_data['<I>/mA'][discharge_data['cycle number'] == n],
                    discharge_data['time/s'][discharge_data['cycle number'] == n]
                )
            )

        discharge = np.abs(np.mean(discharge_list))
        discharge_var = np.abs(np.var(discharge_list))

        ce = np.mean(
            [i / j for i, j in zip(np.abs(discharge_list), charge_list)])
        ce_var = np.var(
            [i / j for i, j in zip(np.abs(discharge_list), charge_list)])

        return charge, discharge, charge_var, discharge_var, ce, ce_var, potential_window

    def create_data_dict(self, file_dict):
        '''
        Creates and then sorts a dictionary of the number of cycles,
        average potentials, average currents, and current variances 
        from the files listed in the file_dict. Entries are grouped
        by scan rate.
        '''

        data_dict = {}

        for key in file_dict:
            data_dict[key] = {'cycles': [],
                              'capacitance': [],
                              'capacitance_var': [],
                              'dis_capacitance': [],
                              'dis_capacitance_var': [],
                              'capacity': [],
                              'capacity_var': [],
                              'dis_capacity': [],
                              'dis_capacity_var': [],
                              'coulombs': [],
                              'coulombs_var': [],
                              'dis_coulombs': [],
                              'dis_coulombs_var': [],
                              'coulomb_eff': [],
                              'c_e_var': []
                              }

        for key in file_dict:

            for element in file_dict[key]:

                with open(element, 'r') as f:
                    data = pd.read_table(f)
                    cycles = max(data['cycle number']) - 2
                    charge, discharge, charge_var, discharge_var, ce, ce_var, potential_window = rate_performance(
                        data=data)

                data_dict[key]['cycles'].append(cycles)

                data_dict[key]['capacitance'].append(charge / potential_window)
                data_dict[key]['capacitance_var'].append(
                    charge_var / potential_window)
                data_dict[key]['dis_capacitance'].append(
                    discharge / potential_window)
                data_dict[key]['dis_capacitance_var'].append(
                    discharge_var / potential_window)

                data_dict[key]['capacity'].append(charge / 3600)
                data_dict[key]['capacity_var'].append(charge_var / 3600)
                data_dict[key]['dis_capacity'].append(discharge / 3600)
                data_dict[key]['dis_capacity_var'].append(discharge_var / 3600)

                data_dict[key]['coulombs'].append(charge / 1000)
                data_dict[key]['coulombs_var'].append(charge_var / 1000)
                data_dict[key]['dis_coulombs'].append(discharge / 1000)
                data_dict[key]['dis_coulombs_var'].append(discharge_var / 1000)

                data_dict[key]['coulomb_eff'].append(ce)
                data_dict[key]['c_e_var'].append(ce_var)

        return data_dict

    def get_weighted_avgs_std(self, data_dict, scalar):
        '''
        Generates dictionary containing the weighted averages of the
        potenials and currents and the current standard deviation for
        each scan rate in the data_dict. Weighted averages use the
        number of cycles recorded at each scan rate
        '''

        scaled_charge_data = np.zeros(len(scan_rates))
        scaled_discharge_data = np.zeros(len(scan_rates))
        scaled_charge_data_std = np.zeros(len(scan_rates))
        scaled_discharge_data_std = np.zeros(len(scan_rates))
        weighted_c_eff = np.zeros(len(scan_rates))
        c_eff_std = np.zeros(len(scan_rates))

        for idx, key in enumerate(data_dict):
            for i in range(len(data_dict[key]['capacitance'])):

                weighted_c_eff[idx] += data_dict[key]['coulomb_eff'][i] * \
                    (data_dict[key]['cycles'][i] /
                     sum(data_dict[key]['cycles']))
                c_eff_std[idx] += data_dict[key]['c_e_var'][i] * \
                    (data_dict[key]['cycles'][i] /
                     sum(data_dict[key]['cycles']))

                if 'Capacitance' in scalar:
                    if 'Specific' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['capacitance'][i] / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['capacitance_var'][i] / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_capacitance'][i] / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_capacitance_var'][i] / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Volumetric' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['capacitance'][i] / (electrode_volume[i] * 1000)) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['capacitance_var'][i] / (
                            electrode_volume[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_capacitance'][i] / (
                            electrode_volume[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_capacitance_var'][i] / (
                            electrode_volume[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Areal' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['capacitance'][i] / (electrode_area[i] * 1000)) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['capacitance_var'][i] / (
                            electrode_area[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_capacitance'][i] / (
                            electrode_area[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_capacitance_var'][i] / (
                            electrode_area[i] * 1000)) * (data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                elif "Capacity" in scalar:
                    if 'Specific' in scalar:
                        scaled_charge_data[idx] += ((data_dict[key]['capacity'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += ((data_dict[key]['capacity_var'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += ((data_dict[key]['dis_capacity'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += ((data_dict[key]['dis_capacity_var'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Volumetric' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['capacity'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['capacity_var'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_capacity'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_capacity_var'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Areal' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['capacity'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['capacity_var'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_capacity'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_capacity_var'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                elif 'Charge' in scalar:
                    if 'Specific' in scalar:
                        scaled_charge_data[idx] += ((data_dict[key]['coulombs'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += ((data_dict[key]['coulombs_var'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += ((data_dict[key]['dis_coulombs'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += ((data_dict[key]['dis_coulombs_var'][i] * 1000) / electrode_mass[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Volumetric' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['coulombs'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['coulombs_var'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_coulombs'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_coulombs_var'][i] / electrode_volume[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

                    elif 'Areal' in scalar:
                        scaled_charge_data[idx] += (data_dict[key]['coulombs'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_charge_data_std[idx] += (data_dict[key]['coulombs_var'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data[idx] += (data_dict[key]['dis_coulombs'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))
                        scaled_discharge_data_std[idx] += (data_dict[key]['dis_coulombs_var'][i] / electrode_area[i]) * (
                            data_dict[key]['cycles'][i] / sum(data_dict[key]['cycles']))

            scaled_charge_data_std[idx] = np.sqrt(scaled_charge_data_std[idx])
            scaled_discharge_data_std[idx] = np.sqrt(
                scaled_discharge_data_std[idx])
            c_eff_std[idx] = np.sqrt(c_eff_std[idx])

        return scaled_charge_data, scaled_discharge_data, scaled_charge_data_std, scaled_discharge_data_std, weighted_c_eff, c_eff_std
