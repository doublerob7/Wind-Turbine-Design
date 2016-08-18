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

data = {}
hist = {}
for month in ['Jan', 'Apr', 'Jul', 'Oct', 'all']:
    data[month] = func.WindDataClass('data/Laramie2005_2015.dat', year='all', month=month)
    hist[month] = list(np.histogram(data[month].speed, bins=25, density=True)[0])

plt.figure()
plt.scatter(list(range(len(hist['Jan']))), hist['Jan'], color='black', marker='x')
plt.scatter(list(range(len(hist['Apr']))), hist['Apr'], color='black', marker='+')
plt.scatter(list(range(len(hist['Jul']))), hist['Jul'], color='black', marker='o')
plt.scatter(list(range(len(hist['Oct']))), hist['Oct'], color='black', marker='x', alpha=0.5)
plt.bar(list(range(len(hist['all']))), hist['all'], color='blue', linestyle='--')
plt.show()




    
''' 3.(b) Using the pdf, determine the mean wind velocity and the wind variance
    for all periods from 3.(a). Provide a table. '''



    
    
''' 3.(c) Overlay a Weibull and Rayleigh distribution on the pdfs and determine
    a best fit for each. '''
    
    
    
    
    
'''  '''