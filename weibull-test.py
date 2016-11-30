import scipy as s
from scipy import stats
# from scipy.stats import exponweib_gen
import functions as f
import matplotlib.pyplot as plt

'''
from http://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm:
x - axis values

Parameters:
Shape - gamma - c
Location - mu
Scale - alpha

If the loc parameter = 0 and the scale parameter = 1, then standard Weibull distribution is given
'''
"""
        Returns
        -------
        shape, loc, scale : tuple of floats
            MLEs for any shape statistics, followed by those for location and
            scale.
"""


data = f.WindDataClass('data/Laramie2005_2015.dat', year='all', month='Jan')

bin_width = 2
data.pdf(bins=range(0, 40, bin_width))

a, shape, loc, scale = stats.exponweib.fit(data.speed, floc=0)

print(a, "shape:", shape, "loc:", loc, "scale:", scale)

x = s.linspace(0, 40, 1000)

plt.hist(data.speed, bins=range(0, 40, bin_width), normed=True)
plt.plot(x, stats.exponweib.pdf(x, a=a, c=shape, scale=scale), 'r--', lw=2, alpha=0.6, label='Exponential Weibull pdf')

plt.show()