# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 20:00:36 2015

ME 4470 - Wind and Ocean Energy
Homework #2

@author: Robert Ressler
"""

''' Setup Litany '''
import numpy as np
import scipy as s
import matplotlib.pyplot as plt

import functions as func

#clear_all()


''' 3.(a) Determine the pdf for the wind velocity for each month, as well as a 
    year. Plot Jan, Apr, July, Oct and a yearly total of all 10 years. '''

jan_data = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month='Jan')
apr_data = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month='Apr')
jul_data = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month='Jul')
oct_data = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month='Oct')
all_data = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month='all')

bin_width = 2
bins = range(0, 40, bin_width)

plt.hist(jan_data.speed, bins, normed=True, histtype='step', color='blue', label='Jan')
plt.hist(apr_data.speed, bins, normed=True, histtype='step', color='red', label='Apr')
plt.hist(jul_data.speed, bins, normed=True, histtype='step', color='orange', label='Jul')
plt.hist(oct_data.speed, bins, normed=True, histtype='step', color='green', label='Oct')
plt.hist(all_data.speed, bins, normed=True, histtype='step', color='black', label='All')

plt.title('Normalized distribution of wind speeds')
plt.grid()
plt.legend()


''' 3.(b) Using the pdf, determine the mean wind velocity and the wind variance
    for all periods from 3.(a). Provide a table. '''

jan_data.pdf(bins=bins)
apr_data.pdf(bins=bins)
jul_data.pdf(bins=bins)
oct_data.pdf(bins=bins)
all_data.pdf(bins=bins)

print("Month |", "Standard Mean |", "PDF Mean |", "Variance")

print('{:^6}'.format(jan_data.month), '{:^15.2f}'.format(jan_data.mean), '{:^10.2f}'.format(jan_data.pdf_mean), '{:^10.2f}'.format(jan_data.variance))
print('{:^6}'.format(apr_data.month), '{:^15.2f}'.format(apr_data.mean), '{:^10.2f}'.format(apr_data.pdf_mean), '{:^10.2f}'.format(apr_data.variance))
print('{:^6}'.format(jul_data.month), '{:^15.2f}'.format(jul_data.mean), '{:^10.2f}'.format(jul_data.pdf_mean), '{:^10.2f}'.format(jul_data.variance))
print('{:^6}'.format(oct_data.month), '{:^15.2f}'.format(oct_data.mean), '{:^10.2f}'.format(oct_data.pdf_mean), '{:^10.2f}'.format(oct_data.variance))
print('{:^6}'.format(all_data.month), '{:^15.2f}'.format(all_data.mean), '{:^10.2f}'.format(all_data.pdf_mean), '{:^10.2f}'.format(all_data.variance))


''' 3.(c) Overlay a Weibull and Rayleigh distribution on the pdfs and determine
    a best fit for each. '''

from scipy.stats import exponweib, rayleigh

x = s.linspace(0, 40, 1000)

plt.figure()
plt.hist(jan_data.speed, bins, normed=True, label='Jan', alpha=0.5, histtype='stepfilled')
a, shape, loc, scale = exponweib.fit(jan_data.speed, floc=0)
plt.plot(x, exponweib.pdf(x, a=a, c=shape, scale=scale), 'b--', label='Exponential Weibull pdf')
loc, scale = rayleigh.fit(jan_data.speed, floc=0)
plt.plot(x, rayleigh.pdf(x, scale=scale), 'r--', label='Rayleigh pdf')
plt.title('Normalized distribution of wind speeds')
plt.grid()

plt.figure()
plt.hist(apr_data.speed, bins, normed=True, label='Apr', alpha=0.6, histtype='stepfilled')
a, shape, loc, scale = exponweib.fit(apr_data.speed, floc=0)
plt.plot(x, exponweib.pdf(x, a=a, c=shape, scale=scale), 'b--', label='Exponential Weibull pdf')
loc, scale = rayleigh.fit(apr_data.speed, floc=0)
plt.plot(x, rayleigh.pdf(x, scale=scale), 'r--', label='Rayleigh pdf')
plt.title('Normalized distribution of wind speeds')
plt.grid()

plt.figure()
plt.hist(jul_data.speed, bins, normed=True, label='Jul', alpha=0.6, histtype='stepfilled')
a, shape, loc, scale = exponweib.fit(jul_data.speed, floc=0)
plt.plot(x, exponweib.pdf(x, a=a, c=shape, scale=scale), 'b--', label='Exponential Weibull pdf')
loc, scale = rayleigh.fit(jul_data.speed, floc=0)
plt.plot(x, rayleigh.pdf(x, scale=scale), 'r--', label='Rayleigh pdf')
plt.title('Normalized distribution of wind speeds')
plt.grid()

plt.figure()
plt.hist(oct_data.speed, bins, normed=True, label='Oct', alpha=0.6, histtype='stepfilled')
a, shape, loc, scale = exponweib.fit(oct_data.speed, floc=0)
plt.plot(x, exponweib.pdf(x, a=a, c=shape, scale=scale), 'b--', label='Exponential Weibull pdf')
loc, scale = rayleigh.fit(jul_data.speed, floc=0)
plt.plot(x, rayleigh.pdf(x, scale=scale), 'r--', label='Rayleigh pdf')
plt.title('Normalized distribution of wind speeds')
plt.grid()

plt.figure()
plt.hist(all_data.speed, bins, normed=True, label='All', alpha=0.6, histtype='stepfilled')
a, shape, loc, scale = exponweib.fit(all_data.speed, floc=0)
plt.plot(x, exponweib.pdf(x, a=a, c=shape, scale=scale), 'b--', label='Exponential Weibull pdf')
loc, scale = rayleigh.fit(jul_data.speed, floc=0)
plt.plot(x, rayleigh.pdf(x, scale=scale), 'r--', label='Rayleigh pdf')
plt.title('Normalized distribution of wind speeds')
plt.grid()

plt.show()
