# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:59 2015

ME 4470 Wind Ocean Energy Homework 5



@author: Robert Ressler
"""

import numpy as np
import matplotlib.pyplot as plot
from Step_3 import TurbineBlade


if __name__ == "__main__":
    """We will again consider the wind turbine blade design we considered in Step 3. Here we will estimate the
    wind turbine power curve by coupling this turbine blade (actually 3 of them) to an induction generator through
    a gearbox. The output goal of this turbine is 425 kW. Torque and power information for the generator is provided
    below. Recall that the turbine was rotating at 30 rpm. However, the input to this generator must be 1800 rpm."""

    from Step_3 import TurbineBlade

    debug = True

    nns = [1.000, 1.005, 1.010, 1.015, 1.020, 1.025, 1.030, 1.035, 1.040]
    T = [0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -2.8, -3.0, -2.8]
    P = [0.0, -100.0, -200.0, -290.0, -375.0, -455.0, -500.0, -500.0, -480.0]


    r = [4, 6, 8, 10, 12, 14, 16, 18, 20]  # m
    R = 21  # m
    chord = [1.5, 1.4, 1.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2]  # m
    theta_p = [26.0, 16.0, 10.0, 5.7, 3.7, 2.5, 2.0, 1.5, 1.0]  # deg
    RPM = 30
    N_blades = 3
    gearratio = 1800 / RPM

    filenames = ("data/DU97W300AoAClCd.csv",) * 3 + ("data/DU91W250AoAClCd.csv",) * 3 + (
                                                                                        "data/DU95W180AoAClCd.csv",) * 3


    """Repeat blade element calculations for 4, 6, 8, 10, and 12 m/s.
    Determine the torque for each of these conditions."""
    ux1 = [4, 6, 8, 10, 12]  # m/s
    turbines = []
    for i, vel in enumerate(ux1):
        turbines.append(TurbineBlade(length=21, min_aero_radius=3, section_radii=r, section_chord=chord,
                         section_twist=theta_p, RPM=RPM, filenames=filenames))
        turbines[i].calculations(ux1=vel, num_blades=N_blades)
        print("{:12}:  {:.2f} kN".format("Thrust", N_blades * turbines[i].thrust * 10 ** -3))
        print("{:12}:  {:.2f} kN m".format("Hub Torque", N_blades * turbines[i].torque * 10 ** -3))
        print("{:12}:  {:.2f} kN m".format("Gen Torque", N_blades * turbines[i].torque * 10 ** -3 / gearratio))
        print("{:12}:  {:.2f}".format("Gear ratio", gearratio), '\n')

    gentorque = []
    for turbine in turbines:
        gentorque.append(turbine.torque / gearratio)

    print(gentorque)

    """Using these torque values, determine the power output from the generator assuming no losses in the bearings
    or gearbox (i.e. power is conserved through the transmission). Present your results as a wind turbine power curve."""

    """If the wind speed got any higher, what would be the danger of continuing to allow the torque to the generator to
    increase? What could be done to ensure the generator still worked effectively at these higher wind speeds?"""

# power = np.zeros(len(ux1))
# torque = np.zeros(len(ux1))
# gentorque = np.zeros(len(ux1))
#
# for i,vel in enumerate(ux1):
#     BEM_data = BEM(r,chord,theta_p,vel,RPM,N_blades,filename)
#     alpha = BEM_data[0]
#     alpha_corr = BEM_data[1]
#     C_N = BEM_data[2]
#     C_T = BEM_data[3]
#     F_N = BEM_data[4]
#     F_T = BEM_data[5]
#     thrust = BEM_data[6]
#     torque[i] = BEM_data[7]
#     power[i] = BEM_data[8]
#     if np.isnan(torque[i]) == True:
#         torque[i] = 0
#     if np.isnan(power[i]) == True:
#         power[i] = 0
#     print 'Hub Torque:  {0:.2f} kN m'.format(torque[i]*10**-3),'Hub Power:  {0:.2f} kW'.format(power[i]*10**-3)
#
#     gentorque[i] = torque[i] / gearratio
#     print 'Generator Torque:  {0:.2f} kN m'.format(-gentorque[i] *10**-3)
#
# #genpower = np.zeros(len(gentorque))
# #for i,tor in enumerate(gentorque):
# #    # genpower[i] is the power produced by the torque fed into the generator
# #    genpower[i] = np.interp(tor,T,P)
#
# plot.figure(1)
# plot.plot(ux1,power*10**-3, 'k--', marker='o', label="Power curve")
# plot.legend(loc='lower right')
# plot.title('Power Curve')
# plot.xlabel('Wind velocity [m/s]')
# plot.ylabel('Power [kW]')
# plot.grid()
    
