# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 22:09:30 2015

@author: robert
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 21:30:21 2015



@author: robert
"""

import numpy as np

# Function to read data file and return parsed data
def ReadNcdcData(filename,yr,mo):
    
    reads = 0    
    badDataCount = 0
    datapoints = 0
    
    with open(filename,'rb') as datafile:
        for i,each in enumerate(datafile):
            reads += 1
            year = each[13:17]
            month = each[17:19]
            
            if year == yr or yr == '':
                if month == mo or mo == '':
                    # First, parse the line and break it into components
                    dttm = each[13:25]
                    print dttm, dttm[0]
                    try:
                        print 'trying to int timemin'
                        timemin = int(dttm[6:8])*1440 + int(dttm[8:10])*60 + int(dttm[10:])
                    except ValueError:
                        print 'failed to create timemin'
                        continue
                    windspeed = each[31:33]
                    winddir = each[26:29]
                    pres = each[107:112]
                    temp = each[85:87]
                    print windspeed, winddir, pres, temp
                    # Next, check for errors or erronious values in the data
                    dataset = [dttm,timemin,windspeed,winddir,pres,temp]
                    print dataset
                    for index, value in enumerate(dataset[2:6],2):
                        try:
                            print datapoints, index, value, value[0]
                            dataset[index] = float(value)
                        except ValueError:
#                            print "Couldn't float("+value+"). Must be a bad data point; skipping this line."
                            badDataCount += 1
                            break
                    else:
                        print dataset
                        try:
                            print 'Appending:'
                            print 'dttm', dataset[0]
                            dttmout.append(dataset[0])
                            print 'timemin',dataset[1]
                            timeminout.append(dataset[1])
                            print 'windspeed',dataset[2]
                            windspeedout.append(dataset[2])
                            print 'winddir', dataset[3]
                            winddirout.append(dataset[3])
                            print 'pres',dataset[4]
                            presout.append(dataset[4])
                            print 'temp',dataset[5]
                            tempout.append(dataset[5])
                            print 'All appends succeeded'
                        except NameError:
                            print 'Excepted'
                            dttmout = [dataset[0]]
                            timeminout = [dataset[1]]
                            windspeedout = [dataset[2]]
                            winddirout = [dataset[3]]
                            presout = [dataset[4]]
                            tempout = [dataset[5]]
                        datapoints += 1
                        
                    continue
                
            
        print "# lines read:",reads,". # of data points:",datapoints,". # data points rejected:",badDataCount,".\n"
    return dttmout,timeminout,windspeedout,winddirout,presout,tempout



[dttm,timemin,windspeed,winddir,pres,temp] = ReadNcdcData('Laramie2005_2015.dat','2005','01')



#yrstart = 0
#mostart = 0
#
#for i,eachdate in enumerate(dttm):
#    year = int(dttm[i][0:4])
#    month = int(dttm[i][4:6]) 
#    if year > int(dttm[yrstart][0:4]):        
#        yrend = i
#        print year-1, yrstart, yrend
#        yrstart = i+1
#        for j,eachyear in enumerate(dttm):
#            if month > int(dttm[mostart][4:6]):
#                moend = j
#                print dttm[mostart], int(dttm[i][0:4]), mostart, moend, np.mean(windspeed[mostart:moend]),
#                mostart = i+1






        