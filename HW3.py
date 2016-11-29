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
import functions as func


class BladeSection:
    """

    """

    def __init__(self, radius, sec_length, chord, twist, filename):
        self.radius = radius
        self.sec_length = sec_length
        self.chord = chord
        self.twist = twist
        self.filename = filename


class Blade:
    """

    """
    def __init__(self, size, sections, section_chord, blade_twist, ):
        self.size = size
        self.sections = sections
        self.chord = section_chord
        self.twist = blade_twist


if __name__ == "__main__":

    """ Given Parameters """

    r = [4,6,8,10,12,14,16,18,20] # m
    chord = [1.5,1.4,1.3,1.2,1.0,0.8,0.6,0.4,0.2] # m
    theta_p = [26,16,10,5.7,3.7,2.5,2.0,1.5,1.0] # deg

    rho = 1.23 # kg/m^3
    ux1 = 12 # m/s

    blade_min_r = 3 # m
    R = 21 # m
    RPM = 30
    omega = (np.float(RPM)/60) * 2 * np.pi # rads/s
    N_blades = 3
    N_sections = 9

    alpha = np.zeros(9)
    alpha_corr = np.zeros(9)
    alpha_corr_deg = np.zeros(9)
    C_N = np.zeros(9)
    C_N_corr = np.zeros(9)
    C_T = np.zeros(9)
    C_T_corr = np.zeros(9)
    F_N = np.zeros(9)
    F_N_corr = np.zeros(9)
    F_T = np.zeros(9)
    F_T_corr = np.zeros(9)

    filename = ['data\DU97W300AoAClCd.csv','data\DU91W250AoAClCd.csv','data\DU95W180AoAClCd.csv']

    # step 1: Guess the induction factors a and a'

    a = 0.1 * np.ones(9)
    a_corr = np.zeros(9)
    a_prime = 0.1 * np.ones(9)
    a_prime_corr = np.zeros(9)

    for i in range(0,8+1):

        tol = 1
        itera = 0

        if i < 3:
            currentfile = filename[0]
        elif i < 6:
            currentfile = filename[1]
        else:
            currentfile = filename[2]

        while tol > 10e-6:

            # step 2: Compute flow angle phi

            tanarg = ( (1-a[i]) * ux1 ) / ( (1+a_prime[i]) * omega * r[i] )
            phi = np.arctan(tanarg)
            phi_deg = np.degrees(phi)

            # step 3: compute local angle of attack alpha

            alpha[i] = phi - np.radians(theta_p[i])
            alpha_deg = np.degrees(alpha)

            # step 4: Determine C_L(alpha) and C_D(alpha) using airfoil properties

            CLCD = func.get_CLCD(currentfile,alpha_deg[i])
            Cl = CLCD[0]
            Cd = CLCD[1]

            # step 5: Determine C_N and C_T from C_L and C_D

            C_N[i] = Cl * np.cos(phi) + Cd * np.sin(phi)
            C_T[i] = Cl * np.sin(phi) - Cd * np.cos(phi)

            # step 6: Calculate a and a'

            s = np.sin(phi)
            c = np.cos(phi)

            sigma = chord[i] * N_blades / (2*np.pi*r[i])

            lasta = a[i]
            a[i] = (1+ (4*(s**2)/(sigma*C_N[i]) ))**(-1)

            lasta_prime = a_prime[i]
            a_prime[i] = ((4*s*c)/(sigma*C_T[i])-1)**-1

            # step 7: Calculate F_N and F_T from C_N and C_T

            w = (ux1 * (1-a[i])) / np.sin(phi)

            F_N[i] = C_N[i] * .5 * rho * chord[i] * (w**2)
            F_T[i] = C_T[i] * .5 * rho * chord[i] * (w**2)

            tol = abs(a[i] - lasta)/lasta
            itera += 1

        print('Sec:',i+1, 'Iter:', itera, ' a: {0:.6}'.format(a[i]),' lasta: {0:.6}'.format(lasta),' alpha: {0:.6}'.format(alpha_deg[i]), ' Tolerance:',tol)

        # Prandtl Tip correction calcs
        f = N_blades/2*(R-r[i])/(r[i]*np.sin(phi))
        F = 2/np.pi * np.arccos(np.exp(-f))

        a_corr[i] = (1+(4*F*s**2)/(sigma*C_N[i]))**-1
        a_prime_corr[i] = ((4*F*s*c)/(sigma*C_T[i])-1)**-1

        phi_corr = np.arctan(((1-a_corr[i])*ux1)/((1+a_prime_corr[i])*omega*r[i]))

        alpha_corr[i] = phi_corr - np.radians(theta_p[i])
        alpha_corr_deg[i] = np.degrees(alpha_corr[i])

        CLCD_corr = func.get_CLCD(currentfile, alpha_corr_deg[i])
        Cl_corr = CLCD_corr[0]
        Cd_corr = CLCD_corr[1]

        C_N_corr[i] = Cl_corr * np.cos(phi_corr) + Cd_corr * np.sin(phi_corr)
        C_T_corr[i] = Cl_corr * np.sin(phi_corr) - Cd_corr * np.cos(phi_corr)

        F_N_corr[i] = C_N_corr[i] * .5 * rho * chord[i] * (w**2)
        F_T_corr[i] = C_T_corr[i] * .5 * rho * chord[i] * (w**2)

        #print '\n',f,F,'\n'

    print('\n')

    """ (a) Plot the local angle of attack, alpha, in degrees as a function of r.
    Plot with and without the Prantl tip correction """

    alpha = np.degrees(alpha)
    alpha_corr = np.degrees(alpha_corr)

    plot.figure(1)
    plot.plot(r,alpha, 'k-', label="AoA without Prandtl correction")
    plot.plot(r,alpha_corr, 'k--', label="AoA with Prandtl correction")
    plot.legend(loc='lower left')
    plot.title('Angle of attack over the length of the blade')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('AoA, alpha [deg]')
    plot.grid()



    """ (b) Plot the local Normal and Tangential force Coefficients as a function
    of r. Plot with and without the Prantl tip correction. """

    plot.figure(2)
    plot.plot(r,C_N_corr, 'k-', label="Normal Coef with Prandtl correction")
    plot.plot(r,C_T_corr, 'k--', label="Tangential Coef with Prandtl correction")
    plot.legend(loc='upper left')
    plot.title('Force Coefficients along blade length')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('Force coefficient')
    plot.grid()

    plot.figure(3)
    plot.plot(r,F_N_corr, 'k-', label="Normal Force with Prandtl correction")
    plot.plot(r,F_T_corr, 'k--', label="Tangential Force with Prandtl correction")
    plot.legend(loc='upper left')
    plot.title('Force along blade length')
    plot.xlabel('Radius: hub to tip midpoints [m]')
    plot.ylabel('Force [N]')
    plot.grid()


    """ (c) Determine total Thrust, Torque and Power experienced by the blades
    using the Prantl correction. """

    # Thrust: the total normal force acting along the entire length of all 3 blades
    sectionlength = r[1] - r[0]
    thrust = 0
    for i, section in enumerate(r):
        thrust += F_N_corr[i] * sectionlength

    thrust *= N_blades
    print('Total Thrust:  {0:.2f} kN'.format(thrust*10**-3))

    # Torque: sum(F_T[i] x r[i]) the sum of each tangential force times it's radius
    torque = 0
    for i, section in enumerate(r):
        torque += F_T_corr[i]*r[i] * sectionlength

    torque *= N_blades
    print('Total Torque:  {0:.2f} kN m'.format(torque*10**-3))

    # Power: torque x angular velocity
    power = torque * omega
    print('Total Power:  {0:.2f} kW'.format(power*10**-3))

    plot.show()
