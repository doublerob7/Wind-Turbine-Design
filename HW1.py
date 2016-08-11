# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 22:48:41 2015

@author: robert
"""

#clear_all()

import calendar

import numpy as np
import matplotlib.pyplot as plt
import functions as func


# Plot the wind speed at 80m from Jan 2014

if 'data' not in locals():
    index = func.index_file('data/Laramie2005_2015.dat')
    data = func.read_wind_data('data/Laramie2005_2015.dat', '2014', 'Jan')

# Convert speeds to 80m
a = 0.19
zref = 10
z = 80

print(max(data.speed))
data.convert_to_hub_height(z, zref, a)
print(max(data.speed))

# plt.figure(1)
# plt.scatter(data.minutes, data.speed, label="Wind Speed at 80m in Jan 2014")
# plt.legend(loc='upper left')
# plt.xlabel('Minutes from beginning of month')
# plt.ylabel('Wind Speed [m/s]')
# plt.grid()
#
# plt.figure(2)
# bins = 40
# plt.hist(data.speed, bins)
# plt.show()


# Plot the average of each month over all years and compare it to the averages from 2014

_2014_data = func.read_wind_data('data/Laramie2005_2015.dat', '2014', 'all')

print(len(_2014_data.speed))
# for month in calendar.month_abbr:
#     if month == '':
#         continue
#     _2014_data.append(func.read_wind_data('data/Laramie2005_2015.dat', '2014', month))

# month_data = {}
# for month in calendar.month_abbr:
#     if month == '':
#         continue
#     for year in range(2005,2015):
#         print('Calling read function with:', year, type(year), month, type(month))
#         try:
#             month_data[str(year)].append(func.read_wind_data('data/Laramie2005_2015.dat', str(year), month))
#         except KeyError:
#             month_data[str(year)] = [func.read_wind_data('data/Laramie2005_2015.dat', str(year), month)]


#plt.figure(2)
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],meanarray, color="blue", label="Monthly mean wind speed 2005-2014")
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],yr2014, color="red", label="Monthly mean wind speed during 2014")
#plt.legend(loc='lower left')
#plt.xlabel('Months')
#plt.ylabel('Mean Wind Speed [m/s]')
