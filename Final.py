# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 10:34:49 2015

ME 4470 Wind and Ocean Energy

Final - Energy Storage

@author: robert
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

	
''' ========================== powerLaw(u,uref,z,zref,alpha) ==================
Normalizes wind speeds from ground height to hub height according to the power 
law.
'''
def powerLaw(uref,z,zref,alpha):
    u = uref * (z/zref) ** alpha
    return u



''' ========================== ReadWindData() =================================
Reads Wind Data from .dat file. Returns relevant data in the following format:

data = [year1, year2, ..., year_n]
year_n = [mo1, mo2, ..., mo12]
mo_n = [dttmarray,timeminarray, windspeedarray, winddirarray, presarray, temparray]

Thus: 
data[0] would be the entire dataset from the first year.
data[0][0] would be the entire dataset from the first month of the first year.
data[0][0][0] would be the first data array from the first month of the first year.
data[0][0][0][0] would be the first entry in the first data array.

EXAMPLE:
data[3][5][2] would be the windspeed array from the 4th year, 6th month. This
could easily be averaged by np.mean(data[3][5][2])
'''
def ReadWindData(filename):
    
    import numpy as np    
    
    datapoints = 0
    badDataCount = 0
    reads = 0
    
    # OPEN FILE
    with open(filename,'rb') as datafile:
        lastdataline = [0,0,0,0,0,0]
        try:
            # LOOP OVER DATA LINES IN DATAFILE
            for dataline in datafile:
                reads +=1
                
                # EXTRACT YEAR, MONTH
                year = dataline[13:17]
                month = dataline[17:19]
                
                # BASIC ERROR HANDLING FOR FIRST LINE
                try:
                    # TRY TO INT THE YEAR, IF THIS FAILS, IT'S NOT A NUMBER
                    year = int(year)
                    month = int(month)
                except ValueError:
                    # IF INT(YEAR) FAILS, CONTINUE TO NEXT LINE IN DATAFILE
                    #print 'int(year) or int(month) failed. Next dataline.'
                    continue
                
                # PARSE THE DATALINE 
                try:
                    dttm = dataline[13:25]
                    timemin = int(dttm[6:8])*1440 + int(dttm[8:10])*60 + int(dttm[10:])
                    windspeed = 0.44704 * np.float(dataline[31:33]) # convert from mph to m/s
                    winddir = np.float(dataline[26:29])
                    if winddir > 360:
                        badDataCount +=1
                        winddir = lastdataline[3]
                
                except:
                    badDataCount +=1
                    windspeed = lastdataline[2]
                    winddir = lastdataline[3]
                        
                
                pres = dataline[107:112]
                temp = dataline[85:87]
                
                # GET THE YEAR POSITION FROM COMPILEDYEARS
                try:
                    # Try to get the position of the current year
                    yearentry = np.where(year==compiledyears)[0][0]
                except:
                    try:
                        # If the value doesn't exist, try to append it
                        compiledyears = np.append(compiledyears,year)
                    except NameError:
                         # Otherwise the array doesn't exist, so create it
                        compiledyears = np.array([year])
                    yearentry = np.where(year==compiledyears)[0][0]
                
                
                if 'data' not in locals():
                    data = [[0]*12]
                
                try:
                    data[yearentry][month-1][0].append(dttm)
                    data[yearentry][month-1][1].append(timemin)
                    data[yearentry][month-1][2].append(windspeed)
                    data[yearentry][month-1][3].append(winddir)
                    data[yearentry][month-1][4].append(pres)
                    data[yearentry][month-1][5].append(temp)
                    datapoints +=1
                except IndexError:
                    data.append([0]*12)
                    data[yearentry][month-1] = [[dttm],[timemin],[windspeed],[winddir],[pres],[temp]]
                    datapoints +=1
                except TypeError:
                    data[yearentry][month-1] = [[dttm],[timemin],[windspeed],[winddir],[pres],[temp]]
                    datapoints +=1
                
                lastdataline = [dttm,timemin,windspeed,winddir,pres,temp]
                
        except:
            print 'Whoa! Something bad happened!'
            print reads, datapoints, dttm
            print compiledyears, year, yearentry
                
    print reads,'Datalines read,',datapoints,'Datapoints kept,',badDataCount,'Datapoints duplicated.'     
    try:
        return data
    except:
        return 0




''' Setup '''
clear_all()
import numpy as np
import matplotlib.pyplot as plot


''' Starting conditions '''
Nturbines = 0
Estored_start = 0


''' Energy consumed per hour '''
Load = [500] * 6 + [1000]*16 + [500]*2
Econs = Load # x 1hr
Econs = Econs[:]


''' Energy produced per hour per turbine'''
powerCurve = [[0,3,4,6,8,10,12,13.6,25.5,25.6],[0.0,0.0,7.3,78.8,190.9,327.6,412.6,425.0,425.0,0.0]]
jandata = ReadWindData('Laramie01_2015.dat')[0][0]
juldata = ReadWindData('Laramie07_2015.dat')[0][6]

for mo in [jandata, juldata]:
    
    # Remove duplicate hour entries
    for i,dttm in enumerate(mo[0][1:len(mo[0])]):
        hour = dttm[8:10]
        if hour == mo[0][i][8:10]:
            try:
                removable.append(i+1)
            except:
                removable = [i+1]
                
    for lists in mo:
        removed = 0
        for each in removable:
            lists.pop(each-removed)
            removed += 1
            
    # Initialize Variables        
    Estored = np.zeros(len(mo[0]))
    Enet = np.zeros(len(mo[2]))
    Econs = np.array(np.tile(Econs,int(mo[0][len(mo[0])-1][6:8])))
    
    # Calculate power at each hour
    for windspeed in mo[2]:
        windspeed = powerLaw(windspeed,80,10,.19)
        
        if windspeed > 25.5:
            windspeed = 0
        
        try:
            power.append(np.interp(windspeed,powerCurve[0],powerCurve[1]))
        except:
            power = [np.interp(windspeed,powerCurve[0],powerCurve[1])]
    
    Eprod = power
    Eprod = np.array(Eprod[:])
    
    # Set the starting point for energy stored
    Estored[0] = Estored_start
    
    
    # (a) Calculate the the number of turbines required to end with a positive value in storage
	
    while Estored[len(Estored)-1] <= 0:
        Nturbines +=1
        Enet = (Eprod * Nturbines) - Econs
        for i,each in enumerate(Enet[0:len(Enet)-1]):
            Estored[i+1] = Estored[i] + each
    
    Egradient = np.gradient(Enet)
    
    Estored_capacity = sum(Enet)
    Estored_start = Estored_capacity
    
    print '\nSolution found! You need', Nturbines,'Turbines.'
    
        
    # (b) Plot power produced and power demanded
    
    plot.figure()
    plot.plot(range(0,len(Eprod)),Eprod*Nturbines, 'k-', label="Power Produced per hour")
    plot.plot(range(0,len(Eprod)),Econs, 'k--', label="Power Demand per hour")
    plot.legend(loc='upper left')
    plot.title('Power Production and Demand')
    plot.xlabel('Hour')
    plot.ylabel('Power [kW]')
    plot.grid()
    
    
    # (c) Determine the minimum capacity of energy storage. Plot energy stored over time
    
    if any(Estored < 0):
        Estored[0] = Estored_start
        for i,each in enumerate(Enet[0:len(Enet)-1]):
            Estored[i+1] = Estored[i] + each
            if Estored[i+1] > Estored_capacity:
                Estored[i+1] = Estored_capacity
    
    Estored_min = min(Estored)
    Estored -= Estored_min
    
    Estored_capacity = Estored_capacity - Estored_min
    Estored_start = Estored_start - Estored_min
    
    plot.figure()
    plot.plot(range(0,len(Eprod)),Enet*10**-3, 'k--', label="Net Power Produced per hour")
    plot.plot(range(0,len(Eprod)),Estored*10**-3, 'k-', label="Energy stored")
    plot.legend(loc='upper left')
    plot.title('Net Power Production and Storage')
    plot.xlabel('Hour')
    plot.ylabel('Power [MW]')
    plot.grid()
    
    print 'Storage Capacity:', int(Estored_capacity)*10**-3 ,'MW, with a starting point of {0:.2f}%'.format(Estored_start/Estored_capacity*100)
    
    
    # (d) Interpret storage capacity in terms of hours on a 500 MW Thermal power plant
    
    ThermalEQhrs = Estored_capacity / (500*10**3)
    print 'That\'s equivalent to {0:.2f} hours on a 500 MW thermal power plant!\n'.format(ThermalEQhrs)
    
    
    
    
    del Eprod, power, Enet, removable
    Estored_start = 0
    Nturbines = 0
    Econs = Load # x 1hr
    Econs = Econs[:]
    