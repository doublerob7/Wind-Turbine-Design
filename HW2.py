# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 20:00:36 2015

ME 4470 - Wind and Ocean Energy
Homework #2

@author: Robert Ressler
"""

''' Setup Litany '''
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt

import functions as func

#clear_all()


''' 3.(a) Determine the pdf for the wind velocity for each month, as well as a 
    year. Plot Jan, Apr, July, Oct and a yearly total of all 10 years. '''

index = func.index_file('data/Laramie2005_2015.dat')
data = {}
for month in ['Jan', 'Apr', 'Jul', 'Oct']:
    data[month] = func.read_wind_data('data/Laramie2005_2015.dat', index, yr='all', mo=month)

bins = range(0, 40+1, 1)

JanWindSpeed = data[0][0][2]
JanWindSpeed = JanWindSpeed[:]
AprWindSpeed = data[0][3][2]
AprWindSpeed = AprWindSpeed[:]
JulWindSpeed = data[0][6][2]
JulWindSpeed = JulWindSpeed[:]
OctWindSpeed = data[0][9][2]
OctWindSpeed = OctWindSpeed[:]
yearWindSpeed = data[0][0][2]
yearWindSpeed = yearWindSpeed[:]

for i in range(0,9+1):
    for speed in data[i][0][2]:
        JanWindSpeed.append(speed)
    for speed in data[i][3][2]:
        AprWindSpeed.append(speed)
    for speed in data[i][6][2]:
        JulWindSpeed.append(speed)
    for speed in data[i][9][2]:
        OctWindSpeed.append(speed)
    for j in range(0,10+1):
        for speed in data[i][j][2]:
            yearWindSpeed.append(speed)
        

Janhistdata = np.histogram(JanWindSpeed, bins, density=True)
Aprhistdata = np.histogram(AprWindSpeed, bins, density=True)
Julhistdata = np.histogram(JulWindSpeed, bins, density=True)
Octhistdata = np.histogram(OctWindSpeed, bins, density=True)
yearhistdata = np.histogram(yearWindSpeed, bins, density=True)

plt.figure(1)
plt.plot(Janhistdata[1][0:len(Janhistdata[1])-1], Janhistdata[0])
plt.figure(2)
plt.plot(Aprhistdata[1][0:len(Janhistdata[1])-1], Aprhistdata[0])
plt.figure(3)
plt.plot(Julhistdata[1][0:len(Janhistdata[1])-1], Julhistdata[0])
plt.figure(4)
plt.plot(Octhistdata[1][0:len(Janhistdata[1])-1], Octhistdata[0])
plt.figure(5)
plt.plot(yearhistdata[1][0:len(Janhistdata[1])-1], yearhistdata[0])
plt.figure(6)
plt.plot(Janhistdata[1][0:len(Janhistdata[1])-1], Janhistdata[0])
plt.plot(Aprhistdata[1][0:len(Janhistdata[1])-1], Aprhistdata[0])
plt.plot(Julhistdata[1][0:len(Janhistdata[1])-1], Julhistdata[0])
plt.plot(Octhistdata[1][0:len(Janhistdata[1])-1], Octhistdata[0])


    
''' 3.(b) Using the pdf, determine the mean wind velocity and the wind variance
    for all periods from 3.(a). Provide a table. '''
    

    
    
''' 3.(c) Overlay a Weibull and Rayleigh distribution on the pdfs and determine
    a best fit for each. '''
    
    
    
    
    
'''  '''