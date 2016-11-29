# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 20:00:36 2015

ME 4470 - Wind and Ocean Energy
Homework #2

@author: Robert Ressler
"""

''' Setup Litany '''
import scipy as s
import matplotlib.pyplot as plt

import functions as func


''' 3.(a) Determine the pdf for the wind velocity for each month, as well as a 
    year. Plot Jan, Apr, July, Oct and a yearly total of all 10 years. '''

bin_width = 2

jan_data = func.WindStats('data/Laramie2005_2015.dat', year='all', month='Jan', bin_width=bin_width)
apr_data = func.WindStats('data/Laramie2005_2015.dat', year='all', month='Apr', bin_width=bin_width)
jul_data = func.WindStats('data/Laramie2005_2015.dat', year='all', month='Jul', bin_width=bin_width)
oct_data = func.WindStats('data/Laramie2005_2015.dat', year='all', month='Oct', bin_width=bin_width)
all_data = func.WindStats('data/Laramie2005_2015.dat', year='all', month='all', bin_width=bin_width)

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

print("Month |", "Standard Mean |", "PDF Mean |", "Variance")

print('{:^6}'.format(jan_data.month), '{:^15.2f}'.format(jan_data.mean), '{:^10.2f}'.format(jan_data.pdf_mean), '{:^10.2f}'.format(jan_data.variance))
print('{:^6}'.format(apr_data.month), '{:^15.2f}'.format(apr_data.mean), '{:^10.2f}'.format(apr_data.pdf_mean), '{:^10.2f}'.format(apr_data.variance))
print('{:^6}'.format(jul_data.month), '{:^15.2f}'.format(jul_data.mean), '{:^10.2f}'.format(jul_data.pdf_mean), '{:^10.2f}'.format(jul_data.variance))
print('{:^6}'.format(oct_data.month), '{:^15.2f}'.format(oct_data.mean), '{:^10.2f}'.format(oct_data.pdf_mean), '{:^10.2f}'.format(oct_data.variance))
print('{:^6}'.format(all_data.month), '{:^15.2f}'.format(all_data.mean), '{:^10.2f}'.format(all_data.pdf_mean), '{:^10.2f}'.format(all_data.variance))


''' 3.(c) Overlay a Weibull and Rayleigh distribution on the pdfs and determine
    a best fit for each. '''

plt.figure()
jan_data.plot_pdf_fit("Jan")
plt.figure()
apr_data.plot_pdf_fit("Apr")
plt.figure()
jul_data.plot_pdf_fit("Jul")
plt.figure()
oct_data.plot_pdf_fit("Oct")
plt.figure()
all_data.plot_pdf_fit("All")

plt.show()
