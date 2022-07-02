import os


class file_reader:

    def __init__(self, rootdir, scan_rates):
        self.rootdir = rootdir
        self.scan_rates = scan_rates

    def create_file_dict(self):
        '''
        Creates a dictionary of file paths that correspond to the cell 
        data that will be processed. Files must be named such that they 
        are pre-sorted. Files are grouped by scan rate.
        '''

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
