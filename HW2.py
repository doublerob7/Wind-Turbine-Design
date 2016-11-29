# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 20:00:36 2015

ME 4470 - Wind and Ocean Energy
Homework #2

@author: Robert Ressler
"""

''' Setup Litany '''
import matplotlib.pyplot as plt
from HW1 import WindDataReader


class WindStats(WindDataReader):
    """ Extends WindDataReader class with analysis methods for wind data """

    def __init__(self, filename, year, month, bin_width, debug=False):
        super().__init__(filename, year, month, debug)
        self._bin_width = bin_width
        self._pdf = None
        self._pdf_mean = None
        self._mean = None
        self._variance = None
        self._weibull_params = None
        self._rayleigh_params = None
        self.bin_centers = []

    @property
    def bin_width(self):
        return self._bin_width

    @bin_width.setter
    def bin_width(self, value):
        if self._bin_width is None:
            self._bin_width = int(value)

    @property
    def variance(self):
        if self._variance is None:
            from scipy import var
            self._variance = var(self.speed)
        return self._variance

    @property
    def mean(self):
        if self._mean is None:
            from scipy import mean
            self._mean = mean(self.speed)
        return self._mean

    @property
    def pdf(self):
        if self._pdf is None:
            from scipy import histogram
            self._pdf = histogram(self.speed, range(0, 50, self.bin_width), density=True)
        return self._pdf

    @property
    def pdf_mean(self):
        if self._pdf_mean is None:
            from scipy import average
            self.calc_bin_centers()
            try:
                self._pdf_mean = average(self.bin_centers, weights=self.pdf[0])
            except:
                print(self.bin_centers, len(self.bin_centers))
                print(self.pdf[0], len(self.pdf[0]))
        return self._pdf_mean

    @property
    def weibull_params(self):
        if self._weibull_params is None:
            from scipy.stats import exponweib
            self._weibull_params = exponweib.fit(self.speed, floc=0)
        return self._weibull_params

    @property
    def rayleigh_params(self):
        if self._rayleigh_params is None:
            from scipy.stats import rayleigh
            self._rayleigh_params = rayleigh.fit(self.speed, floc=0)
        return self._rayleigh_params

    def calc_bin_centers(self):
        center_range = range(0, len(self.pdf[0]))
        for i in center_range:
            next_edge = self.pdf[1][i + 1]
            this_edge = self.pdf[1][i]
            self.bin_centers.append((next_edge - this_edge) / 2 + this_edge)

    def plot_pdf_fit(self, label=None):
        import matplotlib.pyplot as plt
        from scipy.stats import exponweib, rayleigh
        from scipy import linspace, diff

        plt.bar(self.pdf[1][:len(self.pdf[0])], self.pdf[0], width=diff(self.pdf[1]), label=label, alpha=0.5, color='k')

        x = linspace(0, 50, 1000)

        plt.plot(x, exponweib.pdf(x, a=self.weibull_params[0], c=self.weibull_params[1], scale=self.weibull_params[3]), 'b--', label='Exponential Weibull pdf')
        plt.plot(x, rayleigh.pdf(x, scale=self.rayleigh_params[1]), 'r--', label='Rayleigh pdf')

        plt.title('Normalized distribution of wind speeds')
        plt.grid()
        plt.legend()

if __name__ == "__main__":

    ''' 3.(a) Determine the pdf for the wind velocity for each month, as well as a
        year. Plot Jan, Apr, July, Oct and a yearly total of all 10 years. '''

    bin_width = 2

    jan_data = WindStats('data/Laramie2005_2015.dat', year='all', month='Jan', bin_width=bin_width)
    apr_data = WindStats('data/Laramie2005_2015.dat', year='all', month='Apr', bin_width=bin_width)
    jul_data = WindStats('data/Laramie2005_2015.dat', year='all', month='Jul', bin_width=bin_width)
    oct_data = WindStats('data/Laramie2005_2015.dat', year='all', month='Oct', bin_width=bin_width)
    all_data = WindStats('data/Laramie2005_2015.dat', year='all', month='all', bin_width=bin_width)

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
