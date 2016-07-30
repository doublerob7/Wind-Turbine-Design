# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:59 2015

ME 4470 Wind Ocean Energy Homework 5



@author: Robert Ressler
"""

""" Setup """
import numpy as np
import matplotlib.pyplot as plot


nns = [1.000, 1.005, 1.010, 1.015, 1.020, 1.025, 1.030, 1.035, 1.040]
T = [0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -2.8, -3.0, -2.8]
P = [0.0, -100.0, -200.0, -290.0, -375.0, -455.0, -500.0, -500.0, -480.0]

r = [4,6,8,10,12,14,16,18,20] # m
R = 21 # m
chord = [1.5,1.4,1.3,1.2,1.0,0.8,0.6,0.4,0.2] # m
theta_p = [26,16,10,5.7,3.7,2.5,2.0,1.5,1.0] # deg
ux1 = [0,2,4,6,8,10,12] # m/s
RPM = 30
N_blades = 3
filename = ['DU97W300AoAClCd.csv','DU91W250AoAClCd.csv','DU95W180AoAClCd.csv']
gearratio = 1800/RPM

power = np.zeros(len(ux1))
torque = np.zeros(len(ux1))
gentorque = np.zeros(len(ux1))

for i,vel in enumerate(ux1):
    BEM_data = BEM(r,chord,theta_p,vel,RPM,N_blades,filename)
    alpha = BEM_data[0]
    alpha_corr = BEM_data[1]
    C_N = BEM_data[2]
    C_T = BEM_data[3]
    F_N = BEM_data[4]
    F_T = BEM_data[5]
    thrust = BEM_data[6]
    torque[i] = BEM_data[7]
    power[i] = BEM_data[8]
    if np.isnan(torque[i]) == True:
        torque[i] = 0
    if np.isnan(power[i]) == True:
        power[i] = 0
    print 'Hub Torque:  {0:.2f} kN m'.format(torque[i]*10**-3),'Hub Power:  {0:.2f} kW'.format(power[i]*10**-3)
    
    gentorque[i] = torque[i] / gearratio
    print 'Generator Torque:  {0:.2f} kN m'.format(-gentorque[i] *10**-3)

#genpower = np.zeros(len(gentorque))
#for i,tor in enumerate(gentorque):
#    # genpower[i] is the power produced by the torque fed into the generator
#    genpower[i] = np.interp(tor,T,P)
    
plot.figure(1)
plot.plot(ux1,power*10**-3, 'k--', marker='o', label="Power curve")
plot.legend(loc='lower right')
plot.title('Power Curve')
plot.xlabel('Wind velocity [m/s]')
plot.ylabel('Power [kW]')
plot.grid()
    
