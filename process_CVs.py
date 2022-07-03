import numpy as np
import pandas as pd


class process_CVs:
    def __init__(self, mass=None, area=None, volume=None):
        self.mass = mass
        self.area = area
        self.volume = volume

    def cvs_avg_and_variance(self, data):
        """
        Calculates average potential, averager current, and current variance
        for the second to n - 1 cycles for a given scan rate from a BioLogic
        txt file.
        """

        cycle_num = data["cycle number"].unique()

        current = []
        for n in cycle_num[1:-1]:
            current.append(np.asarray(data["<I>/mA"][data["cycle number"] == n]))
        current = pd.DataFrame(current).transpose()
        avg_current = np.mean(current, axis=1)
        current_var = np.var(current, axis=1)

        potential = []
        for n in cycle_num[1:-1]:
            potential.append(np.asarray(data["Ewe/V"][data["cycle number"] == n]))
        potential = pd.DataFrame(potential).transpose()
        avg_potential = potential.mean(axis=1)

        return list(avg_potential), list(avg_current), list(current_var)

    def create_and_sort_data_dict(self, file_dict):
        """
        Creates and then sorts a dictionary of the number of cycles,
        average potentials, average currents, and current variances
        from the files listed in the file_dict. Entries are grouped
        by scan rate.
        """

        data_dict = {}

        for key in file_dict:
            data_dict[key] = {
                "cycles": [],
                "avg_potential": [],
                "avg_current": [],
                "current_var": [],
                "specific_current": [],
                "specific_current_var": [],
                "areal_current": [],
                "areal_current_var": [],
                "volumetric_current": [],
                "volumetric_current_var": [],
            }

        for key in file_dict:

            for element in file_dict[key]:

                with open(element, "r") as f:
                    data = pd.read_table(f)
                    cycles = max(data["cycle number"]) - 2
                    avg_potential, avg_current, current_var = self.cvs_avg_and_variance(
                        data
                    )

                data_dict[key]["cycles"].append(cycles)
                data_dict[key]["avg_potential"].append(avg_potential)
                data_dict[key]["avg_current"].append(avg_current)
                data_dict[key]["current_var"].append(current_var)
                data_dict[key]["specific_current"].append([])
                data_dict[key]["specific_current_var"].append([])
                data_dict[key]["areal_current"].append([])
                data_dict[key]["areal_current_var"].append([])
                data_dict[key]["volumetric_current"].append([])
                data_dict[key]["volumetric_current_var"].append([])

        for key in data_dict:

            for i in range(len(data_dict[key]["avg_potential"])):

                for idx, val in enumerate(data_dict[key]["avg_potential"][i]):

                    if data_dict[key]["avg_potential"][i][0] != max(
                        data_dict[key]["avg_potential"][i]
                    ):

                        pop_potential = data_dict[key]["avg_potential"][i].pop(0)
                        data_dict[key]["avg_potential"][i].append(pop_potential)

                        pop_current = data_dict[key]["avg_current"][i].pop(0)
                        data_dict[key]["avg_current"][i].append(pop_current)

                        pop_current_var = data_dict[key]["current_var"][i].pop(0)
                        data_dict[key]["current_var"][i].append(pop_current_var)

        return data_dict

    def check_and_downsample(self, data_dict):
        """
        Checks number of data points for each entry in the data_dict.
        Downsamples any entry that is longer than the shortet entry.
        """

        lengths = []

        for key in data_dict:
            for idx, val in enumerate(data_dict[key]["avg_potential"]):
                lengths.append(len(val))

        for key in data_dict:
            for idx, val in enumerate(data_dict[key]["avg_potential"]):
                if len(val) != min(lengths):
                    points_to_drop = np.round(
                        np.linspace(
                            0,
                            len(data_dict[key]["avg_potential"][idx]) - 1,
                            len(data_dict[key]["avg_potential"][idx]) - min(lengths),
                        )
                    ).astype(int)

                    potential_update = [
                        data_dict[key]["avg_potential"][idx][i] for i in points_to_drop
                    ]
                    current_update = [
                        data_dict[key]["avg_current"][idx][i] for i in points_to_drop
                    ]
                    var_update = [
                        data_dict[key]["current_var"][idx][i] for i in points_to_drop
                    ]

                    potential_dropped = [
                        item
                        for item in data_dict[key]["avg_potential"][idx]
                        if item not in potential_update
                    ]
                    current_dropped = [
                        item
                        for item in data_dict[key]["avg_current"][idx]
                        if item not in current_update
                    ]
                    var_dropped = [
                        item
                        for item in data_dict[key]["current_var"][idx]
                        if item not in var_update
                    ]

                    data_dict[key]["avg_potential"][idx] = potential_dropped
                    data_dict[key]["avg_current"][idx] = current_dropped
                    data_dict[key]["current_var"][idx] = var_dropped

    def apply_normalization(self, data_dict):
        """
        Normalizes current and current variance according to
        any provided normalization options, i.e., mass, area,
        and/or volume.
        """

        for key in data_dict:
            for idx, val in enumerate(data_dict[key]["avg_potential"]):
                if self.mass != None:
                    data_dict[key]["specific_current"][idx] = [
                        x / self.mass[idx] for x in data_dict[key]["avg_current"][idx]
                    ]
                    data_dict[key]["specific_current_var"][idx] = [
                        x / self.mass[idx] for x in data_dict[key]["current_var"][idx]
                    ]
                if self.area != None:
                    data_dict[key]["areal_current"][idx] = [
                        x / self.area[idx] for x in data_dict[key]["avg_current"][idx]
                    ]
                    data_dict[key]["areal_current_var"][idx] = [
                        x / self.area[idx] for x in data_dict[key]["current_var"][idx]
                    ]
                if self.volume != None:
                    data_dict[key]["volumetric_current"][idx] = [
                        x / self.volume[idx] for x in data_dict[key]["avg_current"][idx]
                    ]
                    data_dict[key]["volumetric_current_var"][idx] = [
                        x / self.volume[idx] for x in data_dict[key]["current_var"][idx]
                    ]

    def get_weighted_avgs_std(self, data_dict):
        """
        Generates dictionary containing the weighted averages of the
        potenials and currents and the current standard deviation for
        each scan rate in the data_dict for each of the normalization
        options that were provided. Weighted averages use the number
        of cycles recorded at each scan rate.
        """

        final_data = {}

        for key in data_dict:
            array_length = len(data_dict[key]["avg_potential"][0])
            potential_window = max(data_dict[key]["avg_potential"][0]) - min(
                data_dict[key]["avg_potential"][0]
            )

            weighted_avg_potential = np.zeros(array_length)
            weighted_avg_current = np.zeros(array_length)
            averaged_current_variance = np.zeros(array_length)

            for idx, val in enumerate(data_dict[key]["avg_potential"]):
                cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                    data_dict[key]["cycles"]
                )
                weighted_avg_potential += np.asarray(val) * cycle_weighting

            for idx, val in enumerate(data_dict[key]["avg_current"]):
                cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                    data_dict[key]["cycles"]
                )
                weighted_avg_current += np.asarray(val) * cycle_weighting

            for idx, val in enumerate(data_dict[key]["current_var"]):
                cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                    data_dict[key]["cycles"]
                )
                averaged_current_variance += np.asarray(val)

            final_data[key] = {
                "w_avg_potential": weighted_avg_potential,
                "w_avg_current": weighted_avg_current,
                "current_std_dev": np.sqrt(averaged_current_variance),
            }

            if self.mass != None:
                weighted_avg_current_density = np.zeros(array_length)
                avg_current_density_variance = np.zeros(array_length)
                weighted_avg_spec_capacitance = np.zeros(array_length)
                avg_spec_capacitance_variance = np.zeros(array_length)
                weighted_avg_spec_capacity = np.zeros(array_length)
                avg_spec_capacity_variance = np.zeros(array_length)
                weighted_avg_spec_charge = np.zeros(array_length)
                avg_spec_charge_variance = np.zeros(array_length)

                for idx, val in enumerate(data_dict[key]["specific_current"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    weighted_avg_current_density += np.asarray(val) * cycle_weighting
                    weighted_avg_spec_capacitance += (
                        np.asarray(val) * cycle_weighting * 1000 / key
                    )
                    weighted_avg_spec_capacity += (
                        np.asarray(val)
                        * cycle_weighting
                        * 1000
                        * potential_window
                        / (3.6 * key)
                    )
                    weighted_avg_spec_charge += (
                        np.asarray(val)
                        * cycle_weighting
                        * 1000
                        * potential_window
                        / key
                    )

                for idx, val in enumerate(data_dict[key]["specific_current_var"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    avg_current_density_variance += np.asarray(val) * cycle_weighting
                    avg_spec_capacitance_variance += np.asarray(val) * cycle_weighting
                    avg_spec_capacity_variance += np.asarray(val) * cycle_weighting
                    avg_spec_charge_variance += np.asarray(val) * cycle_weighting

                final_data[key].update(
                    {
                        "w_avg_current_density": weighted_avg_current_density,
                        "current_density_std_dev": np.sqrt(
                            avg_current_density_variance
                        ),
                        "w_avg_specific_capacitance": weighted_avg_spec_capacitance,
                        "specific_capacitance_std_dev": np.sqrt(
                            avg_spec_capacitance_variance
                        )
                        * 1000
                        / key,
                        "w_avg_specific_capacity": weighted_avg_spec_capacity,
                        "specific_capacity_std_dev": np.sqrt(avg_spec_capacity_variance)
                        * 1000
                        * potential_window
                        / (3.6 * key),
                        "w_avg_specific_charge": weighted_avg_spec_charge,
                        "specific_charge_std_dev": np.sqrt(avg_spec_charge_variance)
                        * 1000
                        * potential_window
                        / key,
                    }
                )

            if self.area != None:
                weighted_avg_areal_current = np.zeros(array_length)
                avg_areal_current_variance = np.zeros(array_length)
                weighted_avg_areal_capacitance = np.zeros(array_length)
                avg_areal_capacitance_variance = np.zeros(array_length)
                weighted_avg_areal_capacity = np.zeros(array_length)
                avg_areal_capacity_variance = np.zeros(array_length)
                weighted_avg_areal_charge = np.zeros(array_length)
                avg_areal_charge_variance = np.zeros(array_length)

                for idx, val in enumerate(data_dict[key]["areal_current"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    weighted_avg_areal_current += np.asarray(val) * cycle_weighting
                    weighted_avg_areal_capacitance += (
                        np.asarray(val) * cycle_weighting / key
                    )
                    weighted_avg_areal_capacity += (
                        np.asarray(val)
                        * cycle_weighting
                        * potential_window
                        / (3.6 * key)
                    )
                    weighted_avg_areal_charge += (
                        np.asarray(val) * cycle_weighting * potential_window / key
                    )

                for idx, val in enumerate(data_dict[key]["areal_current_var"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    avg_areal_current_variance += np.asarray(val) * cycle_weighting
                    avg_areal_capacitance_variance += np.asarray(val) * cycle_weighting
                    avg_areal_capacity_variance += np.asarray(val) * cycle_weighting
                    avg_areal_charge_variance += np.asarray(val) * cycle_weighting

                final_data[key].update(
                    {
                        "w_avg_areal_current": weighted_avg_areal_current,
                        "areal_current_std_dev": np.sqrt(avg_areal_current_variance),
                        "w_avg_areal_capacitance": weighted_avg_areal_capacitance,
                        "areal_capacitance_std_dev": np.sqrt(
                            avg_areal_capacitance_variance
                        )
                        / key,
                        "w_avg_areal_capacity": weighted_avg_areal_capacity,
                        "areal_capacity_std_dev": np.sqrt(avg_areal_capacity_variance)
                        * potential_window
                        / (3.6 * key),
                        "w_avg_areal_charge": weighted_avg_areal_charge,
                        "areal_charge_std_dev": np.sqrt(avg_areal_charge_variance)
                        * potential_window
                        / key,
                    }
                )

            if self.volume != None:
                weighted_avg_volumetric_current = np.zeros(array_length)
                avg_volumetric_current_variance = np.zeros(array_length)
                weighted_avg_volumetric_capacitance = np.zeros(array_length)
                avg_volumetric_capacitance_variance = np.zeros(array_length)
                weighted_avg_volumetric_capacity = np.zeros(array_length)
                avg_volumetric_capacity_variance = np.zeros(array_length)
                weighted_avg_volumetric_charge = np.zeros(array_length)
                avg_volumetric_charge_variance = np.zeros(array_length)

                for idx, val in enumerate(data_dict[key]["volumetric_current"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    weighted_avg_volumetric_current += np.asarray(val) * cycle_weighting
                    weighted_avg_volumetric_capacitance += (
                        np.asarray(val) * cycle_weighting / key
                    )
                    weighted_avg_volumetric_capacity += (
                        np.asarray(val)
                        * cycle_weighting
                        * potential_window
                        / (3.6 * key)
                    )
                    weighted_avg_volumetric_charge += (
                        np.asarray(val) * cycle_weighting * potential_window / key
                    )

                for idx, val in enumerate(data_dict[key]["volumetric_current_var"]):
                    cycle_weighting = data_dict[key]["cycles"][idx] / sum(
                        data_dict[key]["cycles"]
                    )
                    avg_volumetric_current_variance += np.asarray(val) * cycle_weighting
                    avg_volumetric_capacitance_variance += (
                        np.asarray(val) * cycle_weighting
                    )
                    avg_volumetric_capacity_variance += (
                        np.asarray(val) * cycle_weighting
                    )
                    avg_volumetric_charge_variance += np.asarray(val) * cycle_weighting

                final_data[key].update(
                    {
                        "w_avg_volumetric_current": weighted_avg_volumetric_current,
                        "volumetric_current_std_dev": np.sqrt(
                            avg_volumetric_current_variance
                        ),
                        "w_avg_volumetric_capacitance": weighted_avg_volumetric_capacitance,
                        "volumetric_capacitance_std_dev": np.sqrt(
                            avg_volumetric_capacitance_variance
                        )
                        / key,
                        "w_avg_volumetric_capacity": weighted_avg_volumetric_capacity,
                        "volumetric_capacity_std_dev": np.sqrt(
                            avg_volumetric_capacity_variance
                        )
                        * potential_window
                        / (3.6 * key),
                        "w_avg_volumetric_charge": weighted_avg_volumetric_charge,
                        "volumetric_charge_std_dev": np.sqrt(
                            avg_volumetric_charge_variance
                        )
                        * potential_window
                        / key,
                    }
                )

        return final_data
