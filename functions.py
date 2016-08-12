# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:05:45 2015

Utility Functions

@author: robertressler
"""


"""Clears all the variables from the workspace of the spyder application."""
def clear_all():
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]
if __name__ == "__main__":
    clear_all()


class WindDataClass:
    """Contains data and analysis methods for a months worth of wind data.
    """

    def __init__(self, filename, year, month):
        # self.__name__ = name
        self.year = year
        self.month = month
        self.datetime = []
        self.minutes = []
        self.speed = []
        self.direction = []
        self.read_wind_data(filename, year, month)

    def convert_to_hub_height(self, z, zref, a):
        for i, each in enumerate(self.speed):
            self.speed[i] = each * (z / zref) ** a

    def pdf(self, bins='auto'):
        """Calculates the Probability Density Function of the data requested
        """
        import numpy as np
        return np.histogram(self.speed, bins, density=True)

    def read_wind_data(self, filename, yr='', mo=''):
        """Reads Wind Data from .dat 'filename'.
            Returns relevant data by month as a data object of class MonthWindData

            :type filename: str
            :type yr: str
            :type mo: str
            """

        def _suck_data(file, start, end, _reads, _data_points, _badDataCount):
            # LOOP OVER DATA LINES IN DATAFILE
            print(start, end)
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

        import calendar
        import itertools as it

        _months = calendar.month_abbr[1:]
        _data_points = 0
        _badDataCount = 0
        _reads = 0

        # INDEX FILE
        index = index_file(filename)

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
              _badDataCount, 'Datapoints rejected.')



def index_file(filename):
    """Reads a datafile and returns a dict index
    of the starting lines for each year and lists of starting lines for each
     month over each year.
    """

    import calendar

    _reads = 0
    _return_dict = {'keys': []}

    with open(filename, 'r') as _data_file:

        _last_year = ''
        _last_month = ''

        for _data_line in _data_file:

            try:
                _year = str(_data_line[13:17])
                _month = calendar.month_abbr[int(_data_line[17:19])]
            except ValueError:
                continue

            if _year != _last_year:
                _return_dict[_year] = _reads
                _last_year = _year

            if _month != _last_month:
                try:
                    _return_dict[_month].append(_reads)
                except KeyError:
                    _return_dict[_month] = [_reads]
                _last_month = _month

            _reads += 1

    return _return_dict


def read_wind_data(filename, index, yr='all', mo='all'):
    """Reads Wind Data from .dat file.
    Returns relevant data by month as a data object of class MonthWindData

    :type filename: str
    :type index dict
    :type yr: str
    :type mo: str
    """

    import calendar
    import itertools as it

    _months = calendar.month_abbr
    _data_points = 0
    _badDataCount = 0
    _reads = 0
    
    # OPEN FILE
    with open(filename, 'rb') as datafile:

        # SET START AND END POSITIONS FOR DATA BASED ON THE INDEX
        if (yr == 'all') and (mo != 'all'):
            _start_locs = []
            _end_locs = []
            for i, _month in enumerate(_months):
                if _month == '':
                    continue
                _start_locs += [index[_month]]
                _end_locs += [index[_months[i]]]
        elif yr != 'all':
            _start_locs = [index[yr]]
            _end_locs = [index[str(int(yr)+1)]]

        for (start, end) in (_start_locs, _end_locs):
            # LOOP OVER DATA LINES IN DATAFILE
            for _data_line in it.islice(datafile, start, end):
                _reads += 1

                # EXTRACT YEAR, MONTH
                year = str(int(_data_line[13:17]))
                month = _months[int(_data_line[17:19])]

                # PARSE THE DATALINE
                try:
                    winddir = float(_data_line[26:29])
                    if winddir > 360:
                        _badDataCount += 1
                        continue
                    dttm = _data_line[13:25]
                    timemin = int(dttm[6:8])*1440 + int(dttm[8:10])*60 + int(dttm[10:])
                    windspeed = 0.44704 * float(_data_line[31:33])  # convert from mph to m/s

                except ValueError:
                    _badDataCount += 1
                    continue

                # PUT DATA INTO DATA OBJECT
                if '_wind_data' not in locals():
                    _wind_data = MonthWindData(year, month)

                _wind_data.datetime.append(dttm)
                _wind_data.minutes.append(timemin)
                _wind_data.speed.append(windspeed)
                _wind_data.direction.append(winddir)

                _data_points += 1

    print('For ', yr, mo, ': ', _reads, 'Datalines read,', _data_points, 'Datapoints kept,',
          _badDataCount, 'Datapoints rejected.')
    try:
        return _wind_data
    except:
        return None
        

''' ========================== getCLCD(filename,alpha) =================================

Reads C_L and C_D data from 'filename'. Returns interpolated C_L and C_D at alpha.

'''

def getCLCD(filename, alpha):
    
    import numpy as np
    import csv
    #OPEN THE FILE
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        Cl1 = 0
        Cd1 = 0
        alpha1 = 0
        #LOOP OVER DATA LINES
        for _dataline in reader:
            
            #HANDLE THE FIRST LINE THAT ISN'T DATA
            try:
                np.float(_dataline[0])
            except:
                continue
            
            #GRAB THE ALPHA FROM THIS LINE TO COMPARE
            dataalpha = np.float(_dataline[0])
            
            #CHECK ALPHA IN DATA AGAINST ALPHA GIVEN AND STORE Cl AND Cd
            if alpha > dataalpha:
                Cl1 = np.float(_dataline[1])
                Cd1 = np.float(_dataline[2])
                alpha1 = dataalpha
            elif alpha < dataalpha:
                Cl2 = np.float(_dataline[1])
                Cd2 = np.float(_dataline[2])
                
                #CALCULATE SLOPE
                Clm = (Cl2 - Cl1)/(dataalpha - alpha1)
                Cdm = (Cd2 - Cd1)/(dataalpha - alpha1)
                
                #CALCULATE VALUE
                Cl = Cl1 + Clm * (alpha-alpha1)
                Cd = Cd1 + Cdm * (alpha-alpha1)
                
                if Cd1 < 0 or Cd2 < 0:
                    Cd = 0
                
                break
            else:
                Cl = 0
                Cd = 0
                break
    
    try:
        return Cl, Cd
    except:
        return [0, 0]



''' ============ BEM(r[],chord[],theta_p[],ux1,RPM,N_blades,filename['']) =============

BEM() Performs a blade element momentum analysis for given wing parameters. 
Returns alpha, C_N_corr, C_T_corr, F_N_corr, F_T_corr, thrust, torque, power

'''
def BEM(r,chord,theta_p,ux1,RPM,N_blades,filename):
    import numpy as np
    
    R = (r[1]-r[0])/2 + r[len(r)-1]
    rho = 1.23 # kg/m^3
    omega = (np.float(RPM)/60) * 2 * np.pi # rads/s
    
    alpha = np.zeros(9)
    alpha_corr = np.zeros(9)
    alpha_corr_deg = np.zeros(9)
    C_N = np.zeros(9)
    C_N_corr = np.zeros(9)
    C_T = np.zeros(9)
    C_T_corr = np.zeros(9)
    F_N = np.zeros(9)
    F_N_corr = np.zeros(9)
    F_T = np.zeros(9)
    F_T_corr = np.zeros(9) 
    
    # step 1: Guess the induction factors a and a'
        
    a = 0.1 * np.ones(9)
    a_corr = np.zeros(9)
    a_prime = 0.1 * np.ones(9)
    a_prime_corr = np.zeros(9)
    
    for i in range(0,8+1):
    
        tol = 1
        itera = 0
        
        if i < 3:
            currentfile = filename[0]
        elif i < 6:
            currentfile = filename[1]
        else:
            currentfile = filename[2]
        
        while tol > 10e-6:
    
            # step 2: Compute flow angle phi
            
            tanarg = ( (1-a[i])*ux1 ) / ( (1+a_prime[i]) * omega*r[i] )
            phi = np.arctan(tanarg)
            phi_deg = np.degrees(phi)
            
            # step 3: compute local angle of attack alpha
            
            alpha[i] = phi - np.radians(theta_p[i])
            alpha_deg = np.degrees(alpha)
            
            # step 4: Determine C_L(alpha) and C_D(alpha) using airfoil properties
            
            CLCD = getCLCD(currentfile,alpha_deg[i])
            Cl = CLCD[0]
            Cd = CLCD[1]
            
            # step 5: Determine C_N and C_T from C_L and C_D
            
            C_N[i] = Cl * np.cos(phi) + Cd * np.sin(phi)
            C_T[i] = Cl * np.sin(phi) - Cd * np.cos(phi)
            
            # step 6: Calculate a and a'
            
            s = np.sin(phi)
            c = np.cos(phi)
            
            sigma = chord[i] * N_blades / (2*np.pi*r[i])
            
            lasta = a[i]
            a[i] = (1+ (4*(s**2)/(sigma*C_N[i]) ))**(-1)
            
            lasta_prime = a_prime[i]
            a_prime[i] = ((4*s*c)/(sigma*C_T[i])-1)**-1
            
            # step 7: Calculate F_N and F_T from C_N and C_T
            
            w = (ux1 * (1-a[i])) / np.sin(phi) 
            
            F_N[i] = C_N[i] * .5 * rho * chord[i] * (w**2)
            F_T[i] = C_T[i] * .5 * rho * chord[i] * (w**2)
            
            tol = abs(a[i] - lasta)/lasta
            itera += 1
            
        #print 'Sec:',i+1, 'Iter:', itera, ' a: {0:.6}'.format(a[i]),' lasta: {0:.6}'.format(lasta),' alpha: {0:.6}'.format(alpha_deg[i]), ' Tolerance:',tol
            
        # Prandtl Tip correction calcs
        f = N_blades/2*(R-r[i])/(r[i]*np.sin(phi))
        F = 2/np.pi * np.arccos(np.exp(-f))
        
        a_corr[i] = (1+(4*F*s**2)/(sigma*C_N[i]))**-1
        a_prime_corr[i] = ((4*F*s*c)/(sigma*C_T[i])-1)**-1
        
        phi_corr = np.arctan(((1-a_corr[i])*ux1)/((1+a_prime_corr[i])*omega*r[i]))
        
        alpha_corr[i] = phi_corr - np.radians(theta_p[i])
        alpha_corr_deg[i] = np.degrees(alpha_corr[i])
        
        CLCD_corr = getCLCD(currentfile,alpha_corr_deg[i])
        Cl_corr = CLCD_corr[0]
        Cd_corr = CLCD_corr[1]
        
        C_N_corr[i] = Cl_corr * np.cos(phi_corr) + Cd_corr * np.sin(phi_corr)
        C_T_corr[i] = Cl_corr * np.sin(phi_corr) - Cd_corr * np.cos(phi_corr)
        
        F_N_corr[i] = C_N_corr[i] * .5 * rho * chord[i] * (w**2)
        F_T_corr[i] = C_T_corr[i] * .5 * rho * chord[i] * (w**2)
        
    """ Determine total Thrust, Torque and Power experienced by the blades 
    using the Prantl correction. """
    
    #print alpha_deg    
    
    # Thrust: the total normal force acting along the entire length of all 3 blades
    sectionlength = r[1] - r[0]
    thrust = 0
    for i,section in enumerate(r):
        thrust += F_N_corr[i] * sectionlength
            
    thrust *= N_blades
    #print 'Total Thrust:  {0:.2f} kN'.format(thrust*10**-3)
    
    # Torque: sum(F_T[i] x r[i]) the sum of each tangential force times it's radius
    torque = 0
    for i,section in enumerate(r):
        torque += F_T_corr[i]*r[i] * sectionlength
            
    torque *= N_blades
    #print 'Total Torque:  {0:.2f} kN m'.format(torque*10**-3)
    
    # Power: torque x angular velocity
    power = torque * omega
    #print 'Total Power:  {0:.2f} kW'.format(power*10**-3)        
        
    try:
        return np.degrees(alpha),np.degrees(alpha_corr), C_N_corr, C_T_corr, F_N_corr, F_T_corr, thrust, torque, power
    except:
        return [0,0,0,0,0,0,0,0]
        
        
''' ============ BEM(r[],chord[],theta_p[],ux1,RPM,N_blades,filename['']) =============

# good discussion here:  http://stackoverflow.com/questions/4308168/sigmoidal-regression-with-scipy-numpy-python-etc
# curve_fit() example from here: http://permalink.gmane.org/gmane.comp.python.scientific.user/26238
# other sigmoid functions here: http://en.wikipedia.org/wiki/Sigmoid_function

'''

def sigmoid(x, x0, k):
    import numpy as np
    import pylab
    from scipy.optimize import curve_fit
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return y
    



def powerCurve(blade_properties,gen_properties,pitch_theta,vel,RPM):
    """powerCurve takes as input the geometry of a blade (length, chord, twist, blade sections)
    the incoming wind velocity, the hub RPM and the number of blades and returns

    :param blade_properties:
    :param gen_properties:
    :param pitch_theta:
    :param vel:
    :param RPM:
    :return:
    """
    import numpy as np
    
    hubpower = np.zeros(len(vel))
    hubtorque = np.zeros(len(vel))
    gentorque = np.zeros(len(vel))
    
    
    r = blade_properties[0]
    chord = blade_properties[1]
    theta_p = blade_properties[2]
    N_blades = blade_properties[3]
    filename = blade_properties[4]
    theta_p_pitch = np.zeros(len(theta_p))
    
    nns = gen_properties[0]
    T = gen_properties[1]
    P = gen_properties[2]
    gearratio = gen_properties[3]
    
    for i in range(0, len(theta_p)):
        theta_p_pitch[i] = theta_p[i] + pitch_theta
        
    for i, vel in enumerate(vel):
        
        BEM_data = BEM(r, chord, theta_p_pitch, vel, RPM, N_blades, filename)
        
        hubtorque[i] = BEM_data[7]
        hubpower[i] = BEM_data[8]
        
        if np.isnan(hubtorque[i]):
            hubtorque[i] = 0
        if np.isnan(hubpower[i]):
            hubpower[i] = 0
        print('Hub Torque:  {0:.2f} kN m'.format(hubtorque[i]*10**-3), 'Hub Power:  {0:.2f} kW'
              .format(hubpower[i]*10**-3))
        
        if hubpower[i] > 450*10**3:
            for j in range(i, len(hubpower)):
                hubpower[j] = hubpower[i]
            break
        
        gentorque[i] = hubtorque[i] / gearratio
        print('Generator Torque:  {0:.2f} kN m'.format(-gentorque[i] * 10 ** - 3))

    try:
        return [gentorque, hubpower]
        
    except:
        return [0, 0]

