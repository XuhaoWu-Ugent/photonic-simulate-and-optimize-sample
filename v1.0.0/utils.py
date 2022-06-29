from __future__ import print_function
from __future__ import division
import numpy as np
import pandas as pd
import csv


def read_data(filename, column):
    data = pd.read_csv(filename + r'.csv')
    data = data.values
    data_out = data[:, column].transpose()
    return data_out


def write_to_csv(out_name, fieldnames, data):
    if np.shape(data)[1] == len(fieldnames):
        with open(out_name, 'wb') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            arg_str = r'writer.writerow({'
            for idx, fieldname in enumerate(fieldnames):
                arg_str += r'fieldnames[{}]:elem[{}],'.format(idx, idx)
            arg_str += r'})'
            for elem in data:
                exec (arg_str)

    else:
        print("field names has {} columns while data file has {} columns".format(len(fieldnames), np.shape(data)[1]))


def pow2db(pow):
    return 10. * np.log10(pow)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def find_nearest_2d(array, value):
    array = np.asarray(array)
    idx = np.unravel_index((np.abs(array - value)).argmin(), array.shape)
    return array[idx], idx


def find_matching(array, value):
    res = []
    array = np.asarray(array)
    for idx, x in enumerate(np.abs(array - value)):
        if x < 1.0: res.append[idx]
    return res


# For broadband DC coupler L1 and Lt calculation
class Coupler(object):
    """
    Returns a ```Coupler``` object.
    """

    def __init__(self, PhaseCompensationLength):
        self.PhaseCompensationLength = PhaseCompensationLength


class Coupler_DC(Coupler):
    def __init__(self, PhaseCompensationLength, DCLength):
        Coupler.__init__(self, PhaseCompensationLength)
        self.DCLength = DCLength

    def transfer_matrix(self, ne, no, n1, n2, phase1, phase2, wavelength):
        """
        ne: effective index of the even super mode in the DC
        no: effective index of the odd super mode in the DC
        n1: effective index of the fundamental TE/TM mode in the wide waveguide of the phase compensation section
        n2: effective index of the fundamental TE/TM mode in the narrow waveguide of the phase compensation section
        phase1: phase change in the taper joint the DC waveguide to the wide waveguide
        phase2: phase change in the taper joint the DC waveguide to the narrow waveguide
        """
        kappa = np.pi / wavelength * (ne - no)
        k_cross = np.sin(kappa * self.DCLength)
        tau = np.sqrt(1 - k_cross ** 2)
        matrixDC1 = np.array([[tau, -1j * k_cross], [-1j * k_cross, tau]])  # DC 1
        matrixDC2 = matrixDC1  # DC 2
        matrixPC = np.array([[np.exp(-2j * np.pi * n1 * self.PhaseCompensationLength / wavelength), 0],
                             [0, np.exp(-2j * np.pi * n2 * self.PhaseCompensationLength / wavelength)]])
        matrixTaper = np.array([[np.exp(-1j * phase1), 0], [0, np.exp(-1j * phase2)]])
        TM_coupler = np.dot(matrixDC2, np.dot(matrixTaper, np.dot(matrixPC, np.dot(matrixTaper, matrixDC1))))
        return TM_coupler

    def power_cross_single_wavelength(self, ne, no, n1, n2, phase1, phase2, wavelength):
        TM_coupler = self.transfer_matrix(ne, no, n1, n2, phase1, phase2, wavelength)
        return np.abs(TM_coupler[1, 0]) ** 2

    def power_cross(self, polarization):
        """
        read simulated waveguide, DC and taper data from csv file.
        """
        wavelengths, ne, no, n1, n2, phase1, phase2 = read_data(filename=polarization, column=[0, 1, 2, 3, 4, 6, 8])

        TM_coupler_array = []
        for idx, wavelength in enumerate(wavelengths):
            TM_coupler = self.transfer_matrix(ne[idx], no[idx], n1[idx], n2[idx], phase1[idx], phase2[idx], wavelength)
            TM_coupler_array.append(np.abs(TM_coupler[1, 0]) ** 2)
        return np.array(TM_coupler_array)


class Coupler_MMI(Coupler):
    def transfer_matrix(self, n1, n2, phase1, phase2, wavelength):
        """
        n1: effective index of the fundamental TE/TM mode in the wide waveguide of the phase compensation section
        n2: effective index of the fundamental TE/TM mode in the narrow waveguide of the phase compensation section
        phase1: phase change in the taper joint the DC waveguide to the wide waveguide
        phase2: phase change in the taper joint the DC waveguide to the narrow waveguide
        """
        vectorIn = np.array([[np.sqrt(2) / 2], [np.sqrt(2) / 2]])  # DC 1
        vectorOut = vectorIn.transpose()  # DC 2
        matrixPC = np.array([[np.exp(-2j * np.pi * n1 * self.PhaseCompensationLength / wavelength), 0],
                             [0, np.exp(-2j * np.pi * n2 * self.PhaseCompensationLength / wavelength)]])
        matrixTaper = np.array([[np.exp(-1j * phase1), 0], [0, np.exp(-1j * phase2)]])
        TM_coupler = vectorOut.dot(matrixTaper).dot(matrixPC).dot(matrixTaper).dot(vectorIn)
        return TM_coupler

    def power_cross_single_wavelength(self, n1, n2, phase1, phase2, wavelength):
        TM_coupler = self.transfer_matrix(n1, n2, phase1, phase2, wavelength)
        return np.abs(TM_coupler[0, 0]) ** 2

    def power_cross(self, polarization):
        """
        read simulated waveguide, DC and taper data from csv file.
        """
        wavelengths, n1, n2, phase1, phase2 = read_data(filename=polarization, column=[0, 3, 4, 6, 8])

        TM_coupler_array = []
        for idx, wavelength in enumerate(wavelengths):
            TM_coupler = self.transfer_matrix(n1[idx], n2[idx], phase1[idx], phase2[idx], wavelength)
            TM_coupler_array.append(np.abs(TM_coupler[0, 0]) ** 2)
        return np.array(TM_coupler_array)


class Coupler_MMI2X2(Coupler):
    def transfer_matrix(self, n1, n2, phase1, phase2, wavelength):
        """
        n1: effective index of the fundamental TE/TM mode in the wide waveguide of the phase compensation section
        n2: effective index of the fundamental TE/TM mode in the narrow waveguide of the phase compensation section
        phase1: phase change in the taper joint the DC waveguide to the wide waveguide
        phase2: phase change in the taper joint the DC waveguide to the narrow waveguide
        """
        vectorIn = np.array([[np.sqrt(2) / 2], [np.sqrt(2) / 2]])  # MMI1x2 on the input end
        matrixOut = np.sqrt(2) / 2. * np.array([[1, -1j], [-1j, 1]])  # MMI1x2 on the input end
        matrixPC = np.array([[np.exp(-2j * np.pi * n1 * self.PhaseCompensationLength / wavelength), 0],
                             [0, np.exp(-2j * np.pi * n2 * self.PhaseCompensationLength / wavelength)]])
        matrixTaper = np.array([[np.exp(-1j * phase1), 0], [0, np.exp(-1j * phase2)]])
        TM_coupler = matrixOut.dot(matrixTaper).dot(matrixPC).dot(matrixTaper).dot(vectorIn)
        return TM_coupler

    def power_cross_single_wavelength(self, n1, n2, phase1, phase2, wavelength):
        TM_coupler = self.transfer_matrix(n1, n2, phase1, phase2, wavelength)
        return np.abs(TM_coupler[0, 0]) ** 2

    def power_cross(self, polarization):
        """
        read simulated waveguide, DC and taper data from csv file.
        """
        wavelengths, n1, n2, phase1, phase2 = read_data(filename=polarization, column=[0, 3, 4, 6, 8])

        TM_coupler_array = []
        for idx, wavelength in enumerate(wavelengths):
            TM_coupler = self.transfer_matrix(n1[idx], n2[idx], phase1[idx], phase2[idx], wavelength)
            TM_coupler_array.append(np.abs(TM_coupler[0, 0]) ** 2)
        return np.array(TM_coupler_array)


class Coupler_2MMI2X2(Coupler):
    def transfer_matrix(self, n1, n2, phase1, phase2, wavelength):
        """
        n1: effective index of the fundamental TE/TM mode in the wide waveguide of the phase compensation section
        n2: effective index of the fundamental TE/TM mode in the narrow waveguide of the phase compensation section
        phase1: phase change in the taper joint the DC waveguide to the wide waveguide
        phase2: phase change in the taper joint the DC waveguide to the narrow waveguide
        """
        matrixIn = np.sqrt(2) / 2. * np.array([[1, -1j], [-1j, 1]])  # MMI1x2 on the input end
        matrixOut = np.sqrt(2) / 2. * np.array([[1, -1j], [-1j, 1]])  # MMI1x2 on the input end
        matrixPC = np.array([[np.exp(-2j * np.pi * n1 * self.PhaseCompensationLength / wavelength), 0],
                             [0, np.exp(-2j * np.pi * n2 * self.PhaseCompensationLength / wavelength)]])
        matrixTaper = np.array([[np.exp(-1j * phase1), 0], [0, np.exp(-1j * phase2)]])
        TM_coupler = matrixOut.dot(matrixTaper).dot(matrixPC).dot(matrixTaper).dot(matrixIn)
        return TM_coupler

    def power_cross_single_wavelength(self, n1, n2, phase1, phase2, wavelength):
        TM_coupler = self.transfer_matrix(n1, n2, phase1, phase2, wavelength)
        return np.abs(TM_coupler[0, 0]) ** 2

    def power_cross(self, polarization):
        """
        read simulated waveguide, DC and taper data from csv file.
        """
        wavelengths, n1, n2, phase1, phase2 = read_data(filename=polarization, column=[0, 3, 4, 6, 8])

        TM_coupler_array = []
        for idx, wavelength in enumerate(wavelengths):
            TM_coupler = self.transfer_matrix(n1[idx], n2[idx], phase1[idx], phase2[idx], wavelength)
            TM_coupler_array.append(np.abs(TM_coupler[0, 0]) ** 2)
        return np.array(TM_coupler_array)
