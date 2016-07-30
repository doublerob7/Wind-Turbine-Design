# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 22:48:41 2015

@author: robert
"""

%reset

import numpy as np
import matplotlib.pyplot as plt

# Function to read data file and return parsed data
def ReadNcdcData(filename,yr,mo):
      
    badDataCount = 0
    lastbadDataCount = 0
    lastyear = yr
    lastmonth = mo
    
    # OPEN FILE
    with open(filename,'rb') as datafile:
        
        # LOOP OVER DATA LINES IN DATAFILE
        for dataline in datafile:
            
            # EXTRACT YEAR, MONTH
            year = dataline[13:17]
            month = dataline[17:19]
            
            # Control checkpoint
#            print year, [yr], month, [mo]
            
            # BASIC ERROR HANDLING FOR FIRST LINE
            try:
                # TRY TO INT THE YEAR, IF THIS FAILS, IT'S NOT A NUMBER
                int(year)
            except ValueError:
                # IF INT(YEAR) FAILS, CONTINUE TO NEXT LINE IN DATAFILE
                continue
            
            try:
                intmo = int(mo)
                intyr = int(yr)
            except ValueError:
                intmo = 0
                intyr = 0
            
            # Control checkpoint
#            print year, [yr], month, [mo], intmo
            
            # IGNORE DATA WE DON'T WANT
            if year == yr or int(year) == (intyr+1) or yr == '':
#                print '\nPassed year',year, yr, month, [mo], intmo
                if month == mo or int(month) == (intmo+1) or mo == '':
                    pass
#                    print '\nPassed month', month, mo
                else:
#                    print 'Failed month. Continuing...'
                    continue
            else:
                continue
            
            # PARSE THE DATALINE 
            try:
                dttm = dataline[13:25]
                timemin = int(dttm[6:8])*1440 + int(dttm[8:10])*60 + int(dttm[10:])
                windspeed = float(dataline[31:33])
                winddir = float(dataline[26:29])
                pres = float(dataline[107:112])
                temp = float(dataline[85:87])
                if winddir > 360:
                    continue
            except:
                badDataCount +=1
                continue
            
            # Control checkpoint
#            print dttm, timemin, windspeed, winddir, pres, temp
            
            # COLLECT ONE MONTH OF DATA, THEN DO CALCS, STUFF INTO DATA OUTPUT
            try:
                # TRY TO APPEND THE DATALINE TO TEMPDATA
                tempdata.append([dttm, timemin, windspeed, winddir, pres, temp])
                datapoints +=1
#                print "Data appended. Datapoints:", datapoints
            except NameError:
                # IF APPENDING FAILS, THEN IT DOESN'T EXIST. CREATE TEMPDATA
                tempdata = [[dttm, timemin, windspeed, winddir, pres, temp]]
                datapoints = 1
#                print "Appending failed. Created tempdata."
                continue
            
            # Control checkpoint
#            print 'Datapoints:', datapoints
            
            # WHEN COMPILING THE DATA IS DONE, CALCULATE AVERAGES, TRIGGERED BY CHANGE IN MONTH
            lastmonth = tempdata[datapoints-2][0][4:6]
            lastyear = tempdata[datapoints-2][0][0:4]
            if month != lastmonth:
                
#                print 'End of month! Compiling data from:', lastyear, lastmonth
#                print len(tempdata), 'datapoints in month:', lastmonth,'.', badDataCount-lastbadDataCount, 'data points rejected.', '\n'
                lastbadDataCount = badDataCount
                
                # EXTRACT DATA SAMPLE FROM EACH LINE IN TEMPDATA AND THEN AVERAGE IT
                for eachdataline in tempdata:
#                    print eachdataline[eachdataposition]
                    try:
                        dttmarray.append(eachdataline[0])
                        timeminarray.append(eachdataline[1])
                        windspeedarray.append(eachdataline[2])
                        winddirarray.append(eachdataline[3])
                        presarray.append(eachdataline[4])
                        temparray.append(eachdataline[5])
                    except NameError:
                        dttmarray = [eachdataline[0]]
                        timeminarray = [eachdataline[1]]
                        windspeedarray = [eachdataline[2]]
                        winddirarray = [eachdataline[3]]
                        presarray = [eachdataline[4]]
                        temparray = [eachdataline[5]]

                # STUFF EACH DATA ARRAY INTO AN OBJECT                
                averagearray = [windspeedarray,winddirarray,presarray,temparray]
                
                # IF MO IS SPECIFIED, RETURN ALL THE DATA FROM THAT MONTH
                if mo != '':
                    data = tempdata[0:len(tempdata)-1] #[dttmarray,timeminarray,windspeedarray,winddirarray,presarray,temparray]
                
                # IF MO IN UNSPECIFIED, AVERAGE THE DATA TO MAKE ONE LINE PER MONTH
                else:
                    datarow = [lastyear, lastmonth]
                    for j in [0,1,2,3]:
                        datarow.append(np.mean(averagearray[j]))
                
                    # Control checkpoint
    #                print 'dttm:',dttm,'Datarow:',datarow                
                    
                    # PUT THE DATAROW INTO THE DATA OUTPUT OBJECT
                    try:
                        data.append(datarow)
                    except NameError:
                        data = [datarow]
                
                    # CLEAR TEMPDATA BETWEEN EACH MONTH
                    del tempdata, averagearray, datarow
            
            # AFTER COMPILING DATA FOR THE LAST MONTH/YEAR, IGNORE FUTURE ENTRIES
            if yr != '' and year > lastyear:
#                print '\nYear out of bounds:',year, yr
                break
            if mo != '' and month > lastmonth:
#                print '\nMonth out of bounds:',month, mo
                break
    try:
        return data
    except:
        return 0

# 3.a - Plot the wind speed at 80m from Jan 2014

jan = ReadNcdcData('Laramie2005_2015.dat','2014','01')

# Convert speeds to 80m
a = 0.19
zref = 10
z = 80
for entry in jan:
#    print entry[2]
    try:
        speed.append(entry[2])
        timemin.append(entry[1])
    except NameError:
        speed = [entry[2]]
        timemin = [entry[1]]

print np.max(speed)
for i,each in enumerate(speed):
    speed[i] = each*(z/zref)**a
print np.max(speed)

plt.figure(1)
plt.scatter(timemin, speed, label="Wind Speed at 80m in Jan 2014")
plt.legend(loc='upper left')
plt.xlabel('Minutes from beginning of month')
plt.ylabel('Wind Speed [m/s]')

plt.figure(2)
hist = np.histogram(speed, 10, (0,80))
plt.hist(speed)


# 3.b - Plot the average of each month over all years and compare it to the averages from 2014
#data = ReadNcdcData('Laramie2005_2015.dat','','')
#
#print 'mo  mean(velocity)'
#for i,each in enumerate(data):
#    print data[i][0], data[i][1], data[i][2]
#    
#    if each[0] == '2014':
#        try:
#            yr2014.append(each[2])
#        except NameError:
#            yr2014 = [each[2]]
#    
#    if i > 12:
#        i -= 12
#    try:
#        if i == 1:
#            Jan.append(each[2])
#        elif i == 2:
#            Feb.append(each[2])
#        elif i == 3:
#            Mar.append(each[2])
#        elif i == 4:
#            Apr.append(each[2])
#        elif i == 5:
#            May.append(each[2])
#        elif i == 6:
#            Jun.append(each[2])
#        elif i == 7:
#            Jul.append(each[2])
#        elif i == 8:
#            Aug.append(each[2])
#        elif i == 9:
#            Sep.append(each[2])
#        elif i == 10:
#            Oct.append(each[2])
#        elif i == 11:
#            Nov.append(each[2])
#        elif i == 12:
#            Dec.append(each[2])
#    except NameError:
#        if i == 1:
#            Jan = [each[2]]
#        elif i == 2:
#            Feb = [each[2]]
#        elif i == 3:
#            Mar = [each[2]]
#        elif i == 4:
#            Apr = [each[2]]
#        elif i == 5:
#            May = [each[2]]
#        elif i == 6:
#            Jun = [each[2]]
#        elif i == 7:
#            Jul = [each[2]]
#        elif i == 8:
#            Aug = [each[2]]
#        elif i == 9:
#            Sep = [each[2]]
#        elif i == 10:
#            Oct = [each[2]]
#        elif i == 11:
#            Nov = [each[2]]
#        elif i == 12:
#            Dec = [each[2]]
#
#meanarray = [np.mean(Jan),np.mean(Feb),np.mean(Mar),np.mean(Apr),np.mean(May),np.mean(Jun),np.mean(Jul),np.mean(Aug),np.mean(Sep),np.mean(Oct),np.mean(Nov),np.mean(Dec)]
#
#plt.figure(2)
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],meanarray, color="blue", label="Monthly mean wind speed 2005-2014")
#plt.scatter([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],yr2014, color="red", label="Monthly mean wind speed during 2014")
#plt.legend(loc='lower left')
#plt.xlabel('Months')
#plt.ylabel('Mean Wind Speed [m/s]')
