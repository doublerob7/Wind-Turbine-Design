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
    print(index['2014'])
    data = func.read_wind_data('data/Laramie2005_2015.dat', index, '2014', 'Jan')

# Convert speeds to 80m
a = 0.19
zref = 10
z = 80

print(max(data.speed))
data.convert_to_hub_height(z, zref, a)
print(max(data.speed))

plt.figure(1)
plt.scatter(data.minutes, data.speed, label="Wind Speed at 80m in Jan 2014")
plt.legend(loc='upper left')
plt.xlabel('Minutes from beginning of month')
plt.ylabel('Wind Speed [m/s]')
plt.grid()

plt.figure(2)
bins = 40
plt.hist(data.speed, bins)


# Plot the average of each month over all years and compare it to the averages from 2014
_2014_data = {}
month_data = {}
for month in list(calendar.month_abbr):
    if month == '':
        continue
    _2014_data[month] = func.read_wind_data('data/Laramie2005_2015.dat', index, '2014', month)
    for year in range(2005, 2015):
        try:
            _new_month_data = func.read_wind_data('data/Laramie2005_2015.dat', index, str(year), month)
            # print(len(month_data[month]))
            month_data[month] += _new_month_data.speed
            # print('Appended', len(_new_month_data.speed), 'lines. New length:', len(month_data[month]))
        except KeyError:
            month_data[month] = \
                func.read_wind_data('data/Laramie2005_2015.dat', index, str(year), month).speed
            # print('Overwrote speed data')

for month in list(calendar.month_abbr):
    if month == '':
        continue
    try:
        month_data['averages'].append(np.mean(month_data[month]))
        _2014_data['averages'].append(np.mean(_2014_data[month].speed))
    except KeyError:
        month_data['averages'] = [np.mean(month_data[month])]
        _2014_data['averages'] = [np.mean(_2014_data[month].speed)]


plt.figure(3)
plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], month_data['averages'], color="blue", label="Monthly mean wind speed 2005-2014")
plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _2014_data['averages'], color="red", label="Monthly mean wind speed during 2014")
plt.legend(loc='lower left')
plt.xlabel('Months')
plt.ylabel('Mean Wind Speed [m/s]')
plt.show()
