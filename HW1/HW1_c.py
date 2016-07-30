# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 22:48:41 2015

@author: robert
"""

clear_all()

import numpy as np
import matplotlib.pyplot as plt


# 3.a - Plot the wind speed at 80m from Jan 2014

if 'data' not in locals():
    data = ReadWindData('Laramie2005_2015.dat')
jan14 = data[9][0]

# Convert speeds to 80m
a = 0.19
zref = 10
z = 80

speed = jan14[2]
timemin = jan14[1]

print max(speed)
for i,each in enumerate(speed):
    speed[i] = each*(z/zref)**a
print np.max(speed)

plt.figure(1)
plt.scatter(timemin, speed, label="Wind Speed at 80m in Jan 2014")
plt.legend(loc='upper left')
plt.xlabel('Minutes from beginning of month')
plt.ylabel('Wind Speed [m/s]')
plt.grid()

plt.figure(2)
hist = np.histogram(speed, 10, (0,80))
plt.hist(speed)


# 3.b - Plot the average of each month over all years and compare it to the averages from 2014
#data = ReadNcdcData('Laramie2005_2015.dat','','')
#
#print 'mo  mean(velocity)'
#for i,each in enumerate(data):
#    print data[i][0], data[i][1], data[i][2]
#    
#    if each[0] == '2014':
#        try:
#            yr2014.append(each[2])
#        except NameError:
#            yr2014 = [each[2]]
#    
#    if i > 12:
#        i -= 12
#    try:
#        if i == 1:
#            Jan.append(each[2])
#        elif i == 2:
#            Feb.append(each[2])
#        elif i == 3:
#            Mar.append(each[2])
#        elif i == 4:
#            Apr.append(each[2])
#        elif i == 5:
#            May.append(each[2])
#        elif i == 6:
#            Jun.append(each[2])
#        elif i == 7:
#            Jul.append(each[2])
#        elif i == 8:
#            Aug.append(each[2])
#        elif i == 9:
#            Sep.append(each[2])
#        elif i == 10:
#            Oct.append(each[2])
#        elif i == 11:
#            Nov.append(each[2])
#        elif i == 12:
#            Dec.append(each[2])
#    except NameError:
#        if i == 1:
#            Jan = [each[2]]
#        elif i == 2:
#            Feb = [each[2]]
#        elif i == 3:
#            Mar = [each[2]]
#        elif i == 4:
#            Apr = [each[2]]
#        elif i == 5:
#            May = [each[2]]
#        elif i == 6:
#            Jun = [each[2]]
#        elif i == 7:
#            Jul = [each[2]]
#        elif i == 8:
#            Aug = [each[2]]
#        elif i == 9:
#            Sep = [each[2]]
#        elif i == 10:
#            Oct = [each[2]]
#        elif i == 11:
#            Nov = [each[2]]
#        elif i == 12:
#            Dec = [each[2]]
#
#meanarray = [np.mean(Jan),np.mean(Feb),np.mean(Mar),np.mean(Apr),np.mean(May),np.mean(Jun),np.mean(Jul),np.mean(Aug),np.mean(Sep),np.mean(Oct),np.mean(Nov),np.mean(Dec)]
#
#plt.figure(2)
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],meanarray, color="blue", label="Monthly mean wind speed 2005-2014")
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],yr2014, color="red", label="Monthly mean wind speed during 2014")
#plt.legend(loc='lower left')
#plt.xlabel('Months')
#plt.ylabel('Mean Wind Speed [m/s]')
