# -*- coding: utf-8 -*-
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
        
data = ReadWindData('Laramie01_2015.dat')
