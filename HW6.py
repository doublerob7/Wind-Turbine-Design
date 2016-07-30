# -*- coding: utf-8 -*-
"""
ME 4470 Wind and Ocean Energy Homework 6
Created on Wed Dec 02 20:40:55 2015
@author: Robert Ressler
"""

""" Setup """
import numpy as np
import matplotlib.pyplot as plot

''' Definitions '''
nns = [1.000, 1.005, 1.010, 1.015, 1.020, 1.025, 1.030, 1.035, 1.040]
T = [0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -2.8, -3.0, -2.8]
P = [0.0, -100.0, -200.0, -290.0, -375.0, -455.0, -500.0, -500.0, -480.0]

r = [4,6,8,10,12,14,16,18,20] # m
R = 21 # m
chord = [1.5,1.4,1.3,1.2,1.0,0.8,0.6,0.4,0.2] # m
theta_p = [26,16,10,5.7,3.7,2.5,2.0,1.5,1.0] # deg
pitch_theta = [0, 8, 12, 16, 20];
ux1 = range(0,25+1,1) # m/s
RPM = 30
N_blades = 3
filename = ['DU97W300AoAClCd.csv','DU91W250AoAClCd.csv','DU95W180AoAClCd.csv']
gearratio = 1800/30

blade_properties = [r,chord,theta_p,N_blades,filename]
gen_properties = [nns,T,P,gearratio]

''' Preallocation '''
gentorque = np.zeros((len(pitch_theta),len(ux1)))
genpower = np.zeros((len(pitch_theta),len(ux1)))
interpSpeed = np.zeros(len(pitch_theta))


''' Generate a Power Curve for each pitch setting '''
for i,pitch in enumerate(pitch_theta):
    powerData = powerCurve(blade_properties,gen_properties,pitch,ux1,RPM)
    
    gentorque[i] = powerData[0]
    genpower[i] = powerData[1]
    
    # Interpolate for the windspeed corresponding to Rated Power Output
    interpSpeed[i] = np.interp(425*10**3,genpower[i],ux1)

''' Create a new figure '''
plot.figure(1)

''' Plot interpolated intersection points '''
plot.plot(interpSpeed,425*np.ones(len(interpSpeed)), 'ko')

''' Format power curve data to plot only the positive curves up to 450 kW output '''
for powerCurve in genpower:
    for i,power in enumerate(powerCurve):
        
        # trim the negative and 0 values
        if power <= 0:
            continue
        else:
            try:
                # stop compiling data after the first value above 450 kW
                if power == powerCurve[i-1]:
                    break
                posPower.append(power)
                windVel.append(ux1[i])
            except:
                posPower = [power]
                windVel = [ux1[i]]
    # plot the remaining data points
    plot.plot(windVel,[x*10**-3 for x in posPower], 'k--')
    del posPower

''' Format 0 pitch data to extrapolate on calculated interpSpeeds '''
genpowercurve = np.zeros(len(genpower[0]))
genwindvel = np.zeros(len(genpower[0]))
for i,power in enumerate(genpower[0]):
    # collect the power data below 425
    if power < 425*10**3:
        genpowercurve[i] = power
        genwindvel[i] = ux1[i]
    elif power > 425*10**3:
        genpowercurve = genpowercurve[:i]
        genwindvel = genwindvel[:i]

''' Insert the Interpolated wind speeds and powers (rated power) '''
genwindvel = np.insert(genwindvel, len(genwindvel),interpSpeed)
genpowercurve = np.insert(genpowercurve, len(genpowercurve),425*10**3*np.ones(len(interpSpeed)))

''' Plot the composited power curve '''
plot.plot(genwindvel,genpowercurve*10**-3, 'k-', label='Power Curve')
plot.legend(loc='lower right')
plot.title('Power Curve')
plot.xlabel('Wind velocity [m/s]')
plot.ylabel('Power [kW]')
plot.grid()

''' Plot the pitch schedule '''


plot.figure(2)
plot.plot(interpSpeed, pitch_theta, 'ko', label='Pitch')
plot.plot(interpSpeed, pitch_theta, 'k--')
plot.legend(loc='lower right')
plot.title('Pitch Schedule')
plot.xlabel('Wind velocity [m/s]')
plot.ylabel('Pitch angle [deg]')
plot.grid()
