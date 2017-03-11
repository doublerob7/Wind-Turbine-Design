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


class WindDataReader:
    """ Reads a months worth of wind data and parses it. """

    def __init__(self, filename, yr, mo, debug=False):
        self.filename = filename
        self.year = yr
        self.month = mo
        self.datetime = []
        self.minutes = []
        self.speed = []
        self.direction = []
        self.read_wind_data(filename, yr, mo, debug)

    def convert_to_hub_height(self, z=80, zref=10, a=.19):
        for i, each in enumerate(self.speed):
            self.speed[i] = each * (z / zref) ** a

    def index_file(self, filename):
        """Reads a datafile and returns a dict index
        of the starting lines for each year and lists of starting lines for each
         month over each year.
        """

        import calendar

        _reads = 0
        _return_dict = {}

        with open(filename, 'r') as _data_file:

            _last_year = ''
            _last_month = ''

            for _data_line in _data_file:

                _reads += 1

                try:
                    _year = str(_data_line[13:17])
                    _month = calendar.month_abbr[int(_data_line[17:19])]

                    if _year != _last_year:
                        _return_dict[_year] = _reads - 1
                        _last_year = _year

                    if _month != _last_month:
                        try:
                            _return_dict[_month].append(_reads - 1)
                        except KeyError:
                            _return_dict[_month] = [_reads - 1]
                        _last_month = _month

                except ValueError:
                    continue

        return _return_dict

    def read_wind_data(self, filename, yr, mo, debug=False):
        """Reads Wind Data from .dat 'filename'.
            Returns relevant data by month

            :type filename: str
            :type yr: str
            :type mo: str
            :type debug: bool
            """

        def _suck_data(file, start, end, _reads, _data_points, _badDataCount):

            file.seek(0)

            if debug:
                _line = list(it.islice(file, start, start + 1))[0]
                print('For', str(int(_line[13:17])), mo, ', between lines:', start, end)
                print('Parsing', end - start, 'lines. Data is now', len(self.speed), 'lines long')
                _old_data_points = _data_points

                file.seek(0)

            # LOOP OVER DATA LINES IN DATAFILE
            for _data_line in it.islice(file, start, end):
                _reads += 1

                # PARSE THE DATALINE
                try:
                    winddir = float(_data_line[26:29])
                    if winddir > 360:
                        _badDataCount += 1
                        continue
                    dttm = _data_line[13:25]
                    timemin = int(dttm[6:8]) * 1440 + int(dttm[8:10]) * 60 + int(dttm[10:])
                    windspeed = 0.44704 * float(_data_line[31:33])  # convert from mph to m/s

                except ValueError:
                    _badDataCount += 1
                    continue

                # PUT DATA INTO DATA OBJECT
                self.datetime.append(dttm)
                self.minutes.append(timemin)
                self.speed.append(windspeed)
                self.direction.append(winddir)

                _data_points += 1

            if debug:
                print('Parsed', _data_points - _old_data_points, 'new lines. Data is now', len(self.speed),
                  'lines long. Should be:', _data_points, '\n')

            return _reads, _data_points, _badDataCount

        def _next_month(_input):
            """Takes string month input and returns the next month
            :type _input: str
            """
            import calendar
            try:
                return calendar.month_abbr[1:][calendar.month_abbr[1:].index(_input) + 1]
            except IndexError:
                return 'Jan'

        import itertools as it

        _data_points = 0
        _badDataCount = 0
        _reads = 0

        # INDEX FILE
        index = self.index_file(filename)

        # OPEN FILE
        with open(filename, 'rb') as datafile:

            # SET START AND END POSITIONS FOR DATA BASED ON THE INDEX
            if (yr != 'all') and (mo != 'all'):
                _start = index[mo][list(range(2005, 2015)).index(int(yr))]
                _end = index[_next_month(mo)][list(range(2005, 2015)).index(int(yr))]
                (_reads, _data_points, _badDataCount) = \
                    _suck_data(datafile, _start, _end, _reads, _data_points, _badDataCount)

            elif (yr == 'all') and (mo != 'all'):
                _start_locs = index[mo]
                _end_locs = index[_next_month(mo)]

                for (_start, _end) in zip(_start_locs, _end_locs):
                    (_reads, _data_points, _badDataCount) = \
                        _suck_data(datafile, _start, _end, _reads, _data_points, _badDataCount)

            elif (yr != 'all') and (mo == 'all'):
                _start = index[yr]
                _end = index[str(int(yr) + 1)]
                (_reads, _data_points, _badDataCount) = \
                    _suck_data(datafile, _start, _end, _reads, _data_points, _badDataCount)

            else:
                _start = 1
                _end = None
                (_reads, _data_points, _badDataCount) = \
                    _suck_data(datafile, _start, _end, _reads, _data_points, _badDataCount)

        print('For ', yr, mo, ': ', _reads, 'Datalines read,', _data_points, 'Datapoints kept,',
              _badDataCount, 'Datapoints rejected.\n')
        self.convert_to_hub_height()


if __name__ == '__main__':

    # Plot the wind speed at 80m from Jan 2014

    data = WindDataReader('data/Laramie2005_2015.dat', '2014', 'Jan')

    # Convert speeds to 80m
    a = 0.19
    zref = 10
    z = 80

    print(max(data.speed))
    data.convert_to_hub_height(z, zref, a)
    print(max(data.speed))


    plt.figure()
    plt.scatter(data.minutes, data.speed, label="Wind Speed at 80m in Jan 2014")
    plt.legend(loc='upper left')
    plt.xlabel('Minutes from beginning of month')
    plt.ylabel('Wind Speed [m/s]')
    plt.grid()

    # plt.figure(2)
    # bins = 40
    # plt.hist(data.speed, bins)

    # Plot the average of each month over all years and compare it to the averages from 2014
    _2014_data = {}
    month_data = {}
    for month in list(calendar.month_abbr[1:]):
        _2014_data[month] = func.WindDataReader('data/Laramie2005_2015.dat', '2014', month)
        month_data[month] = func.WindDataReader('data/Laramie2005_2015.dat', 'all', month)


    for month in list(calendar.month_abbr[1:]):
        try:
            month_data['averages'].append(np.mean(month_data[month].speed))
            _2014_data['averages'].append(np.mean(_2014_data[month].speed))
        except KeyError:
            month_data['averages'] = [np.mean(month_data[month].speed)]
            _2014_data['averages'] = [np.mean(_2014_data[month].speed)]


    plt.figure()
    plt.plot([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], month_data['averages'], color="blue", label="Monthly mean wind speed 2005-2014")
    plt.plot([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _2014_data['averages'], color="red", label="Monthly mean wind speed during 2014")
    plt.legend(loc='lower left')
    plt.xlabel('Months')
    plt.ylabel('Mean Wind Speed [m/s]')
    plt.grid()
    plt.show()

