# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:59 2015

ME 4470 Wind Ocean Energy Homework 5



@author: Robert Ressler
"""

import numpy as np


class Hub:

    def __init__(self, blade, num_blades):
        self.blade = blade
        self.n_blades = num_blades

    def thrust(self):
        return self.blade.thrust * self.n_blades

    def torque(self):
        return self.blade.torque * self.n_blades

    def power(self):
        return self.blade.power * self.n_blades

    def shaft_speed(self):
        return self.blade.rpm


class GearBox:

    def __init__(self, mechanical_input, ratio, **kwargs):
        self.ratio = ratio
        self.inp = mechanical_input
        if "efficiency" in kwargs:
            self.efficiency = kwargs["efficiency"]
        else:
            self.efficiency = 1

    def torque(self):
        return self.inp.torque() * self.ratio * self.efficiency

    def power(self):
        return self.inp.power() * self.efficiency

    def shaft_speed(self):
        return self.inp.shaft_speed() / self.ratio


class Generator:

    def __init__(self, mechanical_input, slip_list, torque_list, power_list, **kwargs):
        self.slip_list = slip_list
        self.torque_list = tuple(abs(t) for t in torque_list)
        self.power_list = power_list
        self.sync_speed = 1800
        self.rotor = mechanical_input

        if "efficiency" in kwargs:
            self.efficiency = kwargs['efficiency']
        else:
            self.efficiency = 1

    def torque_reaction(self):
        return -np.interp(self.slip(), self.slip_list, self.torque_list)

    def power(self):
        return np.interp(self.slip(), self.slip_list, self.power_list)

    def slip(self):
        torque = self.rotor.torque()
        # print("==================", torque, '\n', self.torque_list, '\n', self.slip_list)
        return np.interp(torque, self.torque_list, self.slip_list)

if __name__ == "__main__":
    """We will again consider the wind turbine blade design we considered in Step 3. Here we will estimate the
    wind turbine power curve by coupling this turbine blade (actually 3 of them) to an induction generator through
    a gearbox. The output goal of this turbine is 425 kW. Torque and power information for the generator is provided
    below. Recall that the turbine was rotating at 30 rpm. However, the input to this generator must be 1800 rpm."""

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

    """Repeat blade element calculations for 4, 6, 8, 10, and 12 m/s.
    Determine the torque for each of these conditions."""
    ux1 = [4, 6, 8, 10, 12, 14, 16, 18, 20]  # m/s
    generated_power = []
    print("Wind Speed: Hub Torque | Generator Rotor Torque")
    for i, vel in enumerate(ux1):
        blade.calculations(ux1=vel, num_blades=3)
        generated_power.append(copy.deepcopy(-generator.power()))
        print("{:6.0f} m/s:{:7.2f} kNm | {:7.2f} kNm".format(vel, hub.torque() * 10 ** -3,
                                                             generator.rotor.torque() * 10 ** -3))

    """Using these torque values, determine the power output from the generator assuming no losses in the bearings
    or gearbox (i.e. power is conserved through the transmission). Present results as a wind turbine power curve."""
    print('\n', "Wind Speed: Gen power")
    for vel, power in zip(ux1, generated_power):
        print("{:6.0f} m/s:{:7.2f} kW".format(vel, power))

    ux1.insert(0, 2)
    generated_power.insert(0, 0)
    ux1.insert(0, 0)
    generated_power.insert(0, 0)

    plot.figure()
    plot.plot(ux1, generated_power, 'k-', marker='o', label="Power curve")
    plot.plot(ux1, [425] * len(ux1), 'k--', label="Rated output")
    plot.legend(loc='lower right')
    plot.title('Power Curve at 0 deg blade pitch')
    plot.xlabel('Wind velocity [m/s]')
    plot.ylabel('Power [kW]')
    plot.grid()

    plot.savefig(filename="output\step5_turbine_power_curve.png", format='png')

    """If the wind speed got any higher, what would be the danger of continuing to allow the torque to the generator to
    increase? What could be done to ensure the generator still worked effectively at these higher wind speeds?"""

    # There's no danger at this blade pitch, as the blade stalls before the hub torque can overwhelm the generator's
    # maximum back-torque.
