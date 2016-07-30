# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:30:17 2015

ME 4470 - Wind and Ocean Energy
Homework #4

@author: Robert Ressler
"""


""" Setup """
import numpy as np
import math
import matplotlib.pyplot as plot

r = [4,6,8,10,12,14,16,18,20] # m
R = 21 # m
chord = [1.5,1.4,1.3,1.2,1.0,0.8,0.6,0.4,0.2] # m
theta_p = [26,16,10,5.7,3.7,2.5,2.0,1.5,1.0] # deg
ux1 = 10 # m/s
RPM = 30
N_blades = 3
filename = ['DU97W300AoAClCd.csv','DU91W250AoAClCd.csv','DU95W180AoAClCd.csv']

t_c = [.3,.3,.3,.25,.25,.25,.18,.18,.18]



# Perform BEM analysis
BEM_data = BEM(r,chord,theta_p,ux1,RPM,N_blades,filename)

alpha = BEM_data[0]
alpha_corr = BEM_data[1]
C_N = BEM_data[2]
C_T = BEM_data[3]
F_N = BEM_data[4]
F_T = BEM_data[5]
thrust = BEM_data[6]
torque = BEM_data[7]
power = BEM_data[8]


# I-beam stats
E = 40*10**9 # Pa
Ix = np.zeros(9)
Iy = np.zeros(9)
A = np.zeros(9)
for i,each in enumerate(r):
    h = .9 * chord[i] * t_c[i]
    b = .15 * chord[i] * t_c[i]
    thickness = .05 * chord[i]
    h1 = h - thickness
    A[i] = b*thickness*2 + h1*thickness
    b_t1 = (b-thickness)/2
    Ix[i] = (1/12)*(b*h**3 - 2*(b_t1*h1**3))
    Iy[i] = ((h * (b**3))/12) - 2*( ((h1*(b_t1**3))/12) + h1 * b_t1 * ( (thickness/2) + (b_t1/2) )**2 )
    T_y = F_N[i] / A
    T_z = F_T[i] / A
    M_y = T_y[i] * r[i]
    M_z = T_z[i] * r[i]
    

""" Plot the shear forces T_y, T_z and the moments M_y and M_z along the blade """
# Shear = V / A, V = F_N, A = from I beam



""" Plot the angular deformations theta_y, theta_z and the deflections u_y and u_z """





