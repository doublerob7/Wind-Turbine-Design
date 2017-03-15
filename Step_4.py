# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:30:17 2015

ME 4470 - Wind and Ocean Energy
Homework #4

@author: Robert Ressler
"""


""" Setup """
import numpy as np
import matplotlib.pyplot as plt

# TODO: implement functions into class
# class Beam:
#     def __init__(self, length, width, height, elastic_mod):
#         self.length = length
#         self.width = (width,)
#         self.height = (height,)
#         self.xsec_area = [width * height for width, height in zip(self.width, self.height)]
#         self.elastic_mod = elastic_mod
#         self.load = None
#         self.support = None
#         self._sections = None
#         self._I = None
#         self._normal = None
#         self._shear = None
#         self._moment = None
#         self._deflection = None
#         self._torsion = None
#         self._twist_angle = None
#
#     def add_load(self, load, loc):
#         try:
#             self.load.append((load, loc))
#         except:
#             self.load = [(load, loc)]
#
#     def add_support(self, type, loc):
#         try:
#             self.support.append((type, loc))
#         except:
#             self.support = [(type, loc)]
#
#     @property
#     def sections(self):
#         if self._sections is None:
#             #calculate it
#             pass
#             self._sections = (((0, 0, 0), (5, 0, 0)), ((5, 0, 0), (10, 0, 0)))
#         return self._sections
#
#     @staticmethod
#     def rect_moment(base, height, d_centroid=0):
#         return (base * (height ** 3)) / 12 + base * height * (d_centroid ** 2)
#
#     @property
#     def I(self):
#         if self._I is None:
#             # Calculate it
#             Ix = []
#             Iy = []
#             for width, height in zip(self.width, self.height):
#                 Ix.append(self.rect_moment(width, height))
#                 Iy.append(self.rect_moment(height, width))
#             self._I = (tuple(Ix), tuple(Iy))
#         return self._I
#
#     @property
#     def shear(self):
#         if self._shear is None:
#             for force in (self.load + self.support):
#                 print(force)
#             self._shear = 0
#         return self._shear


# class IBeam(Beam):
#     def __init__(self, length, width, height, elastic_mod, load, load_loc, t_c):
#         super().__init__(length, width, height, elastic_mod, load, load_loc)
#         self.t_c = t_c
#
#     @property
#     def I(self):
#         if self._I is None:
#             # Calculate it
#             for width, height in zip(self.width, self.height):
#
#                 h = .9 * chord[i] * self.t_c
#                 b = .15 * chord[i] * self.t_c
#                 thickness = .05 * chord[i]
#                 h1 = h - thickness
#                 A[i] = b * thickness * 2 + h1 * thickness
#                 b_t1 = (b - thickness) / 2
#                 self._I_xx = (1 / 12) * (b * h ** 3 - 2 * (b_t1 * h1 ** 3))
#                 self._I_xy = ((h * (b ** 3)) / 12) - 2 * (((h1 * (b_t1 ** 3)) / 12) + h1 * b_t1 * ((thickness / 2) +
#                                                                                                    (b_t1 / 2)) ** 2)
#         return self._I


def trap_integral(x_list, y_list, factor=[1]):
    if factor == [1]:
        factor *= len(y_list)
    values = [0]
    val = 0
    diff_x = list(np.diff(x_list))
    for i, (x_val, y_val, fact) in enumerate(zip(diff_x, y_list, factor)):
        try:
            rect = y_val * x_val
            tri = .5 * (y_list[i+1] - y_val) * x_val
            val += (rect + tri) * fact
            values.append(val)
        except IndexError:
            values.append(0)
        # print(i, x_val, y_val, rect, tri, val, fact)
    return values


def calc_shear(x_list, forces):
    shear_force = [0]
    val = 0
    # print('\n', x, '\n', diff(x), '\n')
    for x_val, force in zip(np.diff(x_list), forces):
        val += force * x_val
        shear_force.append(val)
    return list(reversed(shear_force))


def plot_shear(x_list, shear_list):
    plt.plot(x_list, shear_list, color='b', linewidth=2, label="Shear")
    return plt.plot(x_list, shear_list, color='b', linewidth=2)


def calc_moment(x_list, shear_list):
    _moments = trap_integral(list(val * -1 for val in x_list), list(reversed(shear_list)))
    return list(reversed(_moments))


def plot_moment(x_list, moment_list):
    from matplotlib.pyplot import plot
    plot(x_list, moment_list, color='r', linewidth=2, label="Moment")
    return plot(x_list, moment_list, color='r', linewidth=2)


def moment_inertia(chord, t_c):
    Ix = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    Iy = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    A = [0, 0, 0, 0, 0, 0, 0, 0, 0]

    for i, each in enumerate(r):
        h = .9 * chord[i] * t_c[i]
        b = .15 * chord[i] * t_c[i]
        thickness = .05 * chord[i]
        h1 = h - thickness
        A[i] = b * thickness * 2 + h1 * thickness
        b_t1 = (b - thickness) / 2
        Ix[i] = (1 / 12) * (b * h ** 3 - 2 * (b_t1 * h1 ** 3))
        Iy[i] = ((h * (b ** 3)) / 12) - 2 * (
            ((h1 * (b_t1 ** 3)) / 12) + h1 * b_t1 * ((thickness / 2) + (b_t1 / 2)) ** 2)
    return Ix, Iy


def angular_def(x_list, E, I, moment_list):
    x_rev = list(21 - val for val in reversed(x_list))
    one_over_EI = list(1 / (E * _I) for _I in I)
    _theta = trap_integral(x_rev, moment_list, one_over_EI)
    return _theta


def linear_def(x_list, theta_list):
    x_rev = list(21 - val for val in reversed(x_list))
    _y = trap_integral(x_rev, theta_list)
    return _y


if __name__ == "__main__":

    from Step_3 import TurbineBlade

    debug = True

    r = [4, 6, 8, 10, 12, 14, 16, 18, 20]  # m
    R = 21  # m
    chord = [1.5, 1.4, 1.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2]  # m
    theta_p = [26.0, 16.0, 10.0, 5.7, 3.7, 2.5, 2.0, 1.5, 1.0]  # deg
    ux1 = 10  # m/s
    RPM = 30
    N_blades = 3

    t_c = [.3, .3, .3, .25, .25, .25, .18, .18, .18]

    filenames = ("data/DU97W300AoAClCd.csv",) * 3 + ("data/DU91W250AoAClCd.csv",) * 3 + (
                                                                                        "data/DU95W180AoAClCd.csv",) * 3

    print(filenames)

    blade = TurbineBlade(length=21, min_aero_radius=3, section_radii=r, section_chord=chord,
                         section_twist=theta_p, RPM=RPM, filenames=filenames)

    blade.calculations(ux1=ux1, num_blades=N_blades)

    # I-beam stats
    E = 40 * 10 ** 9  # Pa
    Ix, Iy = moment_inertia(chord, t_c)
    Ix.insert(0, Ix[0])  # duplicate first value to compensate for unloaded section
    Iy.insert(0, Iy[0])

    """ Plot the shear forces T_y, T_z and the moments M_y and M_z along the blade """

    x = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 21]
    load_norm = [force / blade.sec_length for force in blade.sec_norm_force]
    load_norm.reverse()
    load_norm.append(0)

    load_tang = [force / blade.sec_length for force in blade.sec_tang_force]
    load_tang.reverse()
    load_tang.append(0)

    # Shear and moment are calculated from the tip to the hub
    shear_norm = calc_shear(x, load_norm)
    moment_norm = calc_moment(x, shear_norm)

    shear_tang = calc_shear(x, load_tang)
    moment_tang = calc_moment(x, shear_tang)

    if debug:
        print("\n")
        print("{:>15}".format("X:"), len(x), list(reversed(x)))
        print("{:>15}".format("Normal Force:"), len(load_norm), list(reversed(load_norm)))
        print("{:>15}".format("Shear Force:"), len(shear_norm), shear_norm)
        print("{:>15}".format("Moment:"), len(moment_norm), moment_norm)

    plt.figure(1)
    plt.title("Shear and Moment in Normal Direction")
    plot_shear(x, shear_norm)
    plot_moment(x, moment_norm)
    plt.legend()
    plt.grid()
    plt.savefig(filename="output\step4_shear_moment_norm.png", format='png')

    plt.figure(2)
    plt.title("Shear and Moment in Tangential Direction")
    plot_shear(x, shear_tang)
    plot_moment(x, moment_tang)
    plt.legend()
    plt.grid()
    plt.savefig(filename="output\step4_shear_moment_tang.png", format='png')


    """ Plot the angular deformations theta_y, theta_z and the deflections u_y and u_z """
    # dy_dx and dy are calculated from the hub to the tip
    dy_dx = angular_def(x, E, Ix, moment_norm)
    dy = linear_def(x, dy_dx)

    dz_dx = angular_def(x, E, Ix, moment_tang)
    dz = linear_def(x, dz_dx)

    if debug:
        print("{:>15}".format("dy/dx:"), len(dy_dx), dy_dx)
        print("{:>15}".format("y:"), len(dy), dy)
        print("\n")

    plt.figure(3)
    plt.title("Linear Deflections of Blade")
    plt.plot(x, dy, label="y-deflection")
    plt.plot(x, dz, label="z-deflection")
    plt.legend()
    plt.grid()
    plt.savefig(filename="output\step4_blade_deflection.png", format='png')

    # plt.show()
