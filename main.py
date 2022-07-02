from process_CVs import *
from file_reader import *
from rate_performance import *
import matplotlib.pyplot as plt
import logging


def main():

    rootdir = 'C:/Users/Tyler/Desktop/glob_test'
    scan_rates = [0.1, 0.5, 1.0, 2.0, 5.0]

    files = file_reader(rootdir=rootdir, scan_rates=scan_rates)

    file_dict = files.create_file_dict()

    CVs = process_CVs(scan_rates=scan_rates)

    data_dict = CVs.create_and_sort_data_dict(file_dict=file_dict)

    CVs.check_and_downsample(data_dict=data_dict)

    final_data = CVs.get_weighted_avgs_std(data_dict=data_dict)

    plt.plot(final_data[0.1]['w_avg_potential'],
             final_data[0.1]['w_avg_current'])

    plt.show()


if __name__ == "__main__":
    main()
