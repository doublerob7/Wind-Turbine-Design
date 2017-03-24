# -*- coding: utf-8 -*-
"""
ME 4470 Wind and Ocean Energy Homework 6
Created on Wed Dec 02 20:40:55 2015
@author: Robert Ressler
"""

""" Setup """
import numpy as np
from Step_5 import Hub, GearBox, Generator

if __name__ == "__main__":
    """A pitch controlled, constant rotational velocity wind turbine typically operates using one total blade pitch
     theta_p_b_1 (the blade pitch setting at the root) at velocities that produce less than rated power. At higher wind
     velocities, the entire blade is pitched to create increased angles (theta_p) along the blade (this decreases
     the angle of attack) in order to maintain the torque constant such that the generator produces the rated power.
     """

    import matplotlib.pyplot as plot
    from Step_3 import TurbineBlade
    import copy

    debug = True

    # Blade properties
    r = [4, 6, 8, 10, 12, 14, 16, 18, 20]  # m
    R = 21  # m
    chord = [1.5, 1.4, 1.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2]  # m
    theta_p = [26.0, 16.0, 10.0, 5.7, 3.7, 2.5, 2.0, 1.5, 1.0]  # deg
    RPM = 30
    N_blades = 3
    gearratio = 1800 / RPM

    filenames = ("data/DU97W300AoAClCd.csv",) * 3 + ("data/DU91W250AoAClCd.csv",) * 3 + (
                                                                                        "data/DU95W180AoAClCd.csv",) * 3

    # Blade
    blade = TurbineBlade(length=21, min_aero_radius=3, section_radii=r, section_chord=chord,
                         section_twist=theta_p, RPM=RPM, filenames=filenames)

    # Blade Hub
    hub = Hub(blade=blade, num_blades=N_blades)

    # Gearbox
    gearbox = GearBox(mechanical_input=hub, ratio=blade.rpm / 1800)

    # Generator
    nns = [1.000, 1.005, 1.010, 1.015, 1.020, 1.025, 1.030, 1.035]  # dimensionless
    T = [0.0, -500.0, -1000.0, -1500.0, -2000.0, -2500.0, -2800.0, -3000.0]  # N m
    P = [0.0, -100.0, -200.0, -290.0, -375.0, -455.0, -500.0, -500.0]  # kW
    generator = Generator(mechanical_input=gearbox, slip_list=nns, torque_list=T, power_list=P)

    """Now run several velocities for other blade pitches of 8, 12, 16 and 20 degrees that bracket the 425 kW level.
    Be careful not to run velocities that exceed the 425 kW level by too much as your generator will give bad results
    there. Plot these points on the same plot as you plotted your 0 total blade pitch results to produce a result like
    that shown in the figure below."""

    generated_power = {'ux1': [4, 6, 8, 10, 12, 14, 16, 18, 20]}
    power_curve = []
    print("Wind Speed: Hub Torque | Generator Rotor Torque")
    for pitch in [0, 8, 12, 16, 20]:
        generated_power[pitch] = []
        for i, vel in enumerate(generated_power['ux1']):
            blade.calculations(ux1=vel, num_blades=3, pitch=pitch)
            if generator.power() > 0 or np.isnan(generator.power()):
                generated_power[pitch].append(0)
            else:
                generated_power[pitch].append(copy.deepcopy(-generator.power()))
            print("{:6.0f} m/s:  {:7.2f} kW   | {:7.2f} kNm".format(vel, -generator.power(),
                                                                    generator.rotor.torque() * 10 ** -3))

    # power curve
    power_curve = [0, 0] + generated_power[0][:4]
    speeds = [0, 2] + generated_power['ux1'][:4]
    pitch_values = [0] * (len(power_curve) + 1)
    power_curve.append(425)
    speeds.append(np.interp(425, generated_power[0][:6], generated_power['ux1'][:6]))

    for pitch in [8, 12, 16, 20]:
        power_curve.append(425)
        speeds.append(np.interp(425, generated_power[pitch], generated_power['ux1']))
        pitch_values.append(pitch)

    plot.figure()
    title = 'Power Curves with varying blade pitch'
    plot.plot(generated_power['ux1'], [425] * len(generated_power['ux1']), 'k--', label="Rated output")
    plot.plot(generated_power['ux1'][:6], generated_power[0][:6], 'k-.', marker='x', label="0 deg")
    plot.plot(generated_power['ux1'][2:6], generated_power[8][2:6], 'k-.', marker='x', label="8 deg")
    plot.plot(generated_power['ux1'][3:7], generated_power[12][3:7], 'k-.', marker='x', label="12 deg")
    plot.plot(generated_power['ux1'][5:8], generated_power[16][5:8], 'k-.', marker='x', label="16 deg")
    plot.plot(generated_power['ux1'][6:], generated_power[20][6:], 'k-.', marker='x', label="20 deg")
    plot.legend(loc='upper left')
    plot.title(title)
    plot.xlabel('Wind velocity [m/s]')
    plot.ylabel('Power [kW]')
    plot.grid()
    plot.savefig(filename="output\step6_{}.png".format(title), format='png')

    """For each blade setting, estimate the velocity at which the turbine produces rated power. This is the point
    where each total blade pitch curve theta_p_b intersect the rated power. Repeating this for each theta_p_b curve
    yields the circles shown in the figure (plot them on your figure as well). Also generate a new figure that plots
    these blade pitch values against wind speed."""

    plot.figure()
    title = "Turbine System Power Curve"
    plot.plot(speeds, power_curve, 'k-', marker='x', label="System curve")
    plot.plot(speeds, [425] * len(speeds), 'k--', label="Rated output")
    plot.title(title)
    plot.xlabel('Wind velocity [m/s]')
    plot.ylabel('Power [kW]')
    plot.grid()
    plot.savefig(filename="output\step6_{}.png".format(title), format='png')

    plot.figure()
    title = 'Turbine System Pitch control'
    plot.plot(speeds, pitch_values, 'k-', marker='x', label="Pitch value")
    plot.title(title)
    plot.xlabel('Wind velocity [m/s]')
    plot.ylabel('Blade Pitch [deg]')
    plot.grid()
    plot.savefig(filename="output\step6_{}.png".format(title), format='png')

    # plot.show()
