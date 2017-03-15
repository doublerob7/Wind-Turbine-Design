# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 18:49:05 2015

ME 4470 - Wind and Ocean Energy
Homework #3

@author: Robert Ressler
"""

""" Problem 3. Blade Momentum Analysis """

import numpy as np
import matplotlib.pyplot as plot


class TurbineBlade:
    """Class to describe and calculate values for the turbine blade as a whole"""

    def __init__(self, length, min_aero_radius, section_radii, section_twist, section_chord, RPM, filenames):
        from numpy import pi
        self.blade_length = length
        self.min_aero_r = min_aero_radius
        self.sec_radii = section_radii
        self.sec_theta_p = section_twist
        self.sec_chord = section_chord
        self.sec_length = (self.blade_length - self.min_aero_r) / len(self.sec_radii)
        self.rpm = RPM
        self.file_names = filenames
        self.omega = (float(RPM)/60) * 2 * pi  # rads/s
        self.num_sections = len(section_radii)
        self._thrust = None
        self._torque = None
        self._power = None
        self._sec_norm_force = []
        self._sec_tang_force = []
        self._sec_norm_coef = []
        self._sec_tang_coef = []
        self._sec_alpha = []
        self._sec_alpha_corr = []
        self._sec_norm_coef_uncorr = []
        self._sec_tang_coef_uncorr = []

        try:
            assert len(filenames) == len(section_radii)
        except AssertionError:
            print('filenames must be the same length as section_radii')

    @staticmethod
    def get_CLCD(filename, alpha):
        """ Reads C_L and C_D data from 'filename'. Returns interpolated C_L and C_D at alpha.
        :param filename: str
        :param alpha: float
        :return: tuple
        """

        from numpy import float
        import csv
        # OPEN THE FILE
        with open(filename, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            Cl1 = 0
            Cd1 = 0
            alpha1 = 0
            # LOOP OVER DATA LINES
            for _dataline in reader:

                # HANDLE THE FIRST LINE THAT ISN'T DATA
                try:
                    float(_dataline[0])
                except:
                    continue

                # GRAB THE ALPHA FROM THIS LINE TO COMPARE
                dataalpha = float(_dataline[0])

                # CHECK ALPHA IN DATA AGAINST ALPHA GIVEN AND STORE Cl AND Cd
                if alpha > dataalpha:
                    Cl1 = float(_dataline[1])
                    Cd1 = float(_dataline[2])
                    alpha1 = dataalpha
                elif alpha < dataalpha:
                    Cl2 = float(_dataline[1])
                    Cd2 = float(_dataline[2])

                    # CALCULATE SLOPE
                    Clm = (Cl2 - Cl1) / (dataalpha - alpha1)
                    Cdm = (Cd2 - Cd1) / (dataalpha - alpha1)

                    # CALCULATE VALUE
                    Cl = Cl1 + Clm * (alpha - alpha1)
                    Cd = Cd1 + Cdm * (alpha - alpha1)

                    if Cd1 < 0 or Cd2 < 0:
                        Cd = 0

                    break
                else:
                    Cl = 0
                    Cd = 0
                    break

        try:
            return (Cl, Cd)
        except:
            return (0, 0)

    def blade_element_momentum(self, aero_data_file, radius, theta_p, chord, ux1, num_blades, rho, debug=False):
        """Performs iterative BEM analysis on a airfoil section. Returns alpha_deg, alpha_corr_deg, C_N_corr, C_T_corr,
        F_N_corr, F_T_corr.

        :return:
        :type aero_data_file: str
        :type radius: float
        :type theta_p: float
        :type chord: float
        :type ux1: float
        :type num_blades: int
        :type rho: float
        :type debug: bool
        """

        from numpy import arctan, degrees, radians, cos, sin, pi, arccos, exp

        # step 1: Guess the induction factors a and a'
        a = 0.15
        a_prime = 0.15

        tol = 1
        itera = 0

        while tol > 10e-6:

            # step 2: Compute flow angle phi
            phi = arctan(((1 - a) * ux1) / ((1 + a_prime) * self.omega * radius))

            # step 3: compute local angle of attack alpha
            alpha = phi - radians(theta_p)
            alpha_deg = degrees(alpha)

            # step 4: Determine C_L(alpha) and C_D(alpha) using airfoil properties
            Cl, Cd = self.get_CLCD(aero_data_file, alpha_deg)

            # step 5: Determine C_N and C_T from C_L and C_D
            C_N = Cl * cos(phi) + Cd * sin(phi)
            C_T = Cl * sin(phi) - Cd * cos(phi)

            # step 6: Calculate a and a'
            s = sin(phi)
            c = cos(phi)
            sigma = chord * num_blades / (2 * pi * radius)
            lasta = a
            a = (1 + (4 * (s ** 2) / (sigma * C_N))) ** (-1)
            a_prime = ((4 * s * c) / (sigma * C_T) - 1) ** -1

            tol = abs(a - lasta) / lasta
            itera += 1

        if debug:
            print('Sec:', 'Iter:', itera, ' a: {0:.6}'.format(a),' lasta: {0:.6}'.format(lasta),
            ' alpha: {0:.6}'.format(alpha_deg), ' Tolerance:', tol)

        # Prandtl Tip correction calcs
        f = num_blades / 2 * (self.blade_length - radius) / (radius * sin(phi))
        F = 2 / pi * arccos(exp(-f))

        a_corr = (1 + (4 * F * s ** 2) / (sigma * C_N)) ** -1
        a_prime_corr = ((4 * F * s * c) / (sigma * C_T) - 1) ** -1

        phi_corr = arctan(((1 - a_corr) * ux1) / ((1 + a_prime_corr) * self.omega * radius))

        alpha_corr = phi_corr - radians(theta_p)
        alpha_corr_deg = degrees(alpha_corr)

        Cl_corr, Cd_corr = self.get_CLCD(aero_data_file, alpha_corr_deg)

        C_N_corr = Cl_corr * cos(phi_corr) + Cd_corr * sin(phi_corr)
        C_T_corr = Cl_corr * sin(phi_corr) - Cd_corr * cos(phi_corr)

        # step 7: Calculate F_N and F_T from C_N and C_T
        w = (ux1 * (1 - a)) / sin(phi)
        F_N_corr = C_N_corr * .5 * rho * chord * (w ** 2)
        F_T_corr = C_T_corr * .5 * rho * chord * (w ** 2)

        return alpha_deg, alpha_corr_deg, C_N_corr, C_T_corr, F_N_corr, F_T_corr, C_N, C_T

    def calculations(self, ux1=10, num_blades=3, rho=1.23, debug=False):
        for file, radius, theta_P, chord in zip(self.file_names, self.sec_radii, self.sec_theta_p, self.sec_chord):
            a, a_corr, C_N_corr, C_T_corr, F_N, F_T, C_N, C_T = self.blade_element_momentum(file, radius, theta_P, chord, ux1, num_blades, rho, debug)
            self._sec_alpha.append(a)
            self._sec_alpha_corr.append(a_corr)
            self._sec_norm_coef.append(C_N_corr)
            self._sec_tang_coef.append(C_T_corr)
            self._sec_norm_force.append(F_N)
            self._sec_tang_force.append(F_T)
            self._sec_norm_coef_uncorr.append(C_N)
            self._sec_tang_coef_uncorr.append(C_T)
        return self

    @property
    def thrust(self):
        if self._thrust is None:
            self._thrust = sum((force * self.sec_length for force in self.sec_norm_force))
        return self._thrust

    @property
    def torque(self):
        if self._torque is None:
            self._torque = sum((force * radius * self.sec_length for force, radius in zip(self.sec_tang_force, self.sec_radii)))
        return self._torque

    @property
    def power(self):
        if self._power is None:
            self._power = self.torque * self.omega
        return self._power

    @property
    def sec_norm_force(self):
        if self._sec_norm_force is []:
            self.calculations()
        return self._sec_norm_force

    @property
    def sec_tang_force(self):
        if self._sec_tang_force is []:
            self.calculations()
        return self._sec_tang_force

    @property
    def sec_norm_coef(self):
        if self._sec_norm_coef is []:
            self.calculations()
        return self._sec_norm_coef

    @property
    def sec_tang_coef(self):
        if self._sec_tang_coef is []:
            self.calculations()
        return self._sec_tang_coef

    @property
    def sec_alpha(self):
        if self._sec_alpha is []:
            self.calculations()
        return self._sec_alpha

    @property
    def sec_alpha_corr(self):
        if self._sec_alpha_corr is []:
            self.calculations()
        return self._sec_alpha_corr


if __name__ == '__main__':

    """ Given Parameters """
    r = [4, 6, 8, 10, 12, 14, 16, 18, 20]  # m
    chord = [1.5, 1.4, 1.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2]  # m
    theta_p = [26, 16, 10, 5.7, 3.7, 2.5, 2.0, 1.5, 1.0]  # deg

    rho = 1.23  # kg/m^3
    ux1 = 10  # m/s

    blade_min_r = 3  # m
    R = 21  # m
    RPM = 30
    omega = (np.float(RPM)/60) * 2 * np.pi  # rads/s
    N_blades = 3
    N_sections = 9

    filenames = ['data\DU97W300AoAClCd.csv'] * 3 + ['data\DU91W250AoAClCd.csv'] * 3 + ['data\DU95W180AoAClCd.csv'] * 3
    blade = TurbineBlade(length=R, min_aero_radius=blade_min_r, section_radii=r, section_twist=theta_p,
                         section_chord=chord, RPM=RPM, filenames=filenames)

    blade.calculations(ux1=ux1, num_blades=N_blades, rho=rho, debug=False)


    """ (a) Plot the local angle of attack, alpha, in degrees as a function of r.
    Plot with and without the Prantl tip correction """

    plot.figure(1)
    plot.plot(blade.sec_radii, blade.sec_alpha, 'k-', label="AoA without Prandtl correction")
    plot.plot(blade.sec_radii, blade.sec_alpha_corr, 'k--', label="AoA with Prandtl correction")
    plot.legend(loc='upper left')
    plot.title('Angle of attack over the length of the blade')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('AoA, alpha [deg]')
    plot.grid()
    plot.savefig(filename="output\step3_AoA.png", format='png')

    """ (b) Plot the local Normal and Tangential force Coefficients as a function
    of r. Plot with and without the Prantl tip correction. """

    plot.figure(num=2, figsize=(7, 8))
    plot.subplot(2, 1, 1)
    plot.plot(blade.sec_radii, blade._sec_norm_coef_uncorr, 'k--', label="Normal Coef without Prandtl correction")
    plot.plot(blade.sec_radii, blade.sec_norm_coef, 'k-', label="Normal Coef with Prandtl correction")
    plot.legend(loc='lower right')
    plot.title('Force Coefficients along blade length')
    plot.ylabel('Force coefficient')
    plot.grid()
    plot.subplot(2, 1, 2)
    plot.plot(blade.sec_radii, blade._sec_tang_coef_uncorr, 'k--', label="Tangential Coef without Prandtl correction")
    plot.plot(blade.sec_radii, blade.sec_tang_coef, 'k-', label="Tangential Coef with Prandtl correction")
    plot.legend(loc='upper right')
    plot.title('Force Coefficients along blade length')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('Force coefficient')
    plot.grid()
    plot.savefig(filename="output\step3_force_coeff.png", format='png')

    plot.figure(3)
    plot.plot(blade.sec_radii, blade.sec_norm_force, 'k-', label="Normal Force")
    plot.plot(blade.sec_radii, blade.sec_tang_force, 'k--', label="Tangential Force")
    plot.legend(loc='upper left')
    plot.title('Force along blade length')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('Force [N]')
    plot.grid()
    plot.savefig(filename="output\step3_blade_loads.png", format='png')


    """ (c) Determine total Thrust, Torque and Power experienced by the blades
    using the Prantl correction. """
    print('Total Thrust:  {0:.2f} kN'.format(N_blades * blade.thrust * 10 ** -3))
    print('Total Torque:  {0:.2f} kN m'.format(N_blades * blade.torque * 10 ** -3))
    print('Total Power:  {0:.2f} kW'.format(N_blades * blade.power * 10 ** -3))

    with open("output\step3_hub_parameters.txt", "w") as file:
        # Thrust: the total normal force acting along the entire length of all 3 blades
        print('Total Thrust:  {0:.2f} kN'.format(N_blades * blade.thrust * 10 ** -3), file=file)

        # Torque: sum(F_T[i] x r[i]) the sum of each tangential force times it's radius
        print('Total Torque:  {0:.2f} kN m'.format(N_blades * blade.torque * 10 ** -3), file=file)

        # Power: torque x angular velocity
        print('Total Power:  {0:.2f} kW'.format(N_blades * blade.power * 10 ** -3), file=file)

    # plot.show()
