import datetime
import csv
import numpy as np
import sys

from movementStrategies import getDatesBetweenRange
## Read in files ##
def createDataMatrix(filename, resolution):
    """Transforms text file into usable matrix. Text files downloaded from https://lter.limnology.wisc.edu/

    Args:
        filename: the text file 
    Returns: 
        dataMatrix: list of lists
    """
    dataMatrix = []
    errorFile = filename + "error.txt"
    f = open(errorFile, "w")
    
    with open(filename, "r") as ins:
        f.write(next(ins))
        for line in ins: 
            row = line.strip().split(',')
            savedData = []
            try:
                if resolution == 'hourly':
                    #These are the headers:
                    #sampledate,year4,month,daynum,hour,depth,wtemp,flag_wtemp
                    sampleDate = formatDateHourly(row[0]) #sampledate
                    sampleHour = int(float(row[4])/100) # hour, corrects for military time eg. 0100 -> 1
                    
                    savedData.append(sampleDate.replace(hour = sampleHour)) #sets to correct time
                    savedData.append(float(row[5])) # depth
                    savedData.append(float(row[6])) #wtemp
                    savedData.append(row[7]) #flag_wtemp
                    dataMatrix.append(savedData) 
               
                elif resolution == 'daily':
                    #These are the headers:
                    #sampledate,year4,month,daynum,depth,wtemp,flag_wtemp
                    sampleDate = dateString(row[0]) #sampledate, this sets fake hour & min

                    savedData.append(sampleDate)
                    savedData.append(float(row[4])) # depth
                    savedData.append(float(row[5])) #wtemp
                    savedData.append(row[6]) #flag_wtemp
                    dataMatrix.append(savedData) 

                elif resolution == 'hires':
                    #These are the headers:
                    #sampledate,year4,month,daynum,sample_time,data_freq,depth,wtemp,flag_wtemp
                    sampleDate = formatDateHourly(row[0]) #sampleDate
                    sampleTime = row[4] #sample_time
                    sampleHour = sampleTime[:2]
                    sampleMin = sampleTime[2:]

                    savedData.append(sampleDate.replace(hour = sampleHour, min = sampleMin)) #sets to correct time
                    savedData.append(float(row[6])) #depth
                    savedData.append(float(row[7])) #wtemp
                    savedData.append(row[8]) #flag_wtemp
                    savedData.append(row[5]) #data_freq

            except:
                f.write(','.join(row)+'\n')
               
    f.close()
    return dataMatrix

def createSunriseData(filename, year):
    """Transforms sunrise data from a plain text file to a useable matrix. 
    The textfiles used were manually currated from https://aa.usno.navy.mil/index.php 
    using Wausau, Wisconsin as the location. 

    Args:
        filename: the text file 
        year: the year of interest

    Returns: 
        sunriseData: a matrix containing all of the sunrise and sunset data for that year
    """
    sunriseData = [year]
    errorFile = filename + "error.txt"
    f = open(errorFile, "w")
    
    with open(filename, "r") as ins:
        f.write(next(ins))
        f.write(next(ins))
        f.write(next(ins))
        for line in ins:
            try:
                row = []
                row.append(int(line[:2]))
                line = line[4:]
                line = [line[i:i+11] for i in range(0, len(line), 11)] #this is the best way to account for uneven whitespace when missing values
                
                for l in line:
                    if l.isspace():
                        row.append([])
                    else:
                        row.append([int(i) for i in l.strip().split()])

                sunriseData.append(row)
            except:
               f.write(','.join(row)+'\n') 
    return sunriseData      

## convert and transform data ##

def createDictionary(dataMatrix):
    """Transforms dataMatrix into a nested dictionary.

    Args:
        dataMatrix: list of lists
    Returns: 
        dataMatrixDict: nested dictionary 
    """
    dataMatrixDict = {}
    for element in dataMatrix:
        dateKey = element[0]
        if dateKey not in dataMatrixDict:
            dataMatrixDict[dateKey] = {}
        if element[1] not in dataMatrixDict[dateKey]:
            dataMatrixDict[dateKey][element[1]] = [element[2], element[3]]

    return dataMatrixDict

def groupAllDates(dataMatrixDict):
    """ Group all dates found in dataMatrixDict by their month/day.

        Args: 
            dataMatrixDict: the dictionary of depths and temperatures

        Returns:
            dateDictionary: grouped dictionary with each key as a Month-Day string
    """
    dateDictionary = {}
    for date in sorted(dataMatrixDict.keys()):
        dateKey = date.strftime('%m-%d %H:%M:%S')
        if dateKey not in dateDictionary:
            dateDictionary[dateKey] = {}
        for depth in dataMatrixDict[date]:
            if depth not in dateDictionary[dateKey]:
                dateDictionary[dateKey][depth] = []
            dateDictionary[dateKey][depth].append(dataMatrixDict[date][depth][0])
   
    return dateDictionary

def averageOverYears(dateDictionary):
    """Finds the average temperature at each point over all years.

    Args:
        dateDictionary: grouped dictionary with each key as a Month-Day string

    Returns:
        dateDictionary: grouped dictionary with lists of real temps replaced with single mean temp
    """
    for date in dateDictionary.keys():
        for depth in dateDictionary[date]:
            dateDictionary[date][depth] = np.mean(dateDictionary[date][depth])

    return dateDictionary

def fillGapsInData(dataMatrixDict, dateDictionary, singleYear):
    """Fills each gap in time with an averaged value across all years.
    
    Args: 
        dataMatrixDict: the dictionary of depths and temperatures
        dateDictionary: grouped dictionary with lists of real temps replaced with single mean temp
        singleYear: the year of interest
    
    Returns:
        dataMatrixDict: final dictionary with gaps filled
    """
    missingHours = 0
    missingDepths = 0
    actualDays = sorted(dataMatrixDict.keys())

    allDepths = set()

    first, last = dateString(singleYear, True)

    allDays = getDatesBetweenRange(first, last, True)
    for date in dataMatrixDict.keys():
        if date in allDays:
            allDepths.update(set(dataMatrixDict[date].keys()))

    counter = 0
    for date in allDays:
        counter += 1
        if counter%500 == 0:
            print '%d of %d'%(counter,len(allDays))
        dateKey = date.strftime('%m-%d %H:%M:%S')
        if dateKey not in dateDictionary:
            continue #because we don't have an average for it
        if date not in actualDays:
            dataMatrixDict[date] = {}
            for depth in allDepths:
                if depth in dateDictionary[dateKey]:
                    missingHours += 1
                    dataMatrixDict[date][depth] = [dateDictionary[dateKey][depth]]
        else:
            for depth in allDepths:
                if depth not in dataMatrixDict[date] and depth in dateDictionary[dateKey]:
                    missingDepths += 1
                    dataMatrixDict[date][depth] = [dateDictionary[dateKey][depth]]

    print "missing Hours: %d \nmissing Depths: %d"%(missingHours, missingDepths)
    return dataMatrixDict

def extendDataMatrix(dataMatrixDict, startDate, endDate, resolution=0.1,verbose=False):
    """Handles data interpolation. 

    Args:
        dataMatrixDict: nested dictionary 
        startDate: the start date
        endDate: the end date
    Returns: 
        dataMatrixDict: interpolated nested dictionary
    """
    myDate = startDate

    while myDate <= endDate:
        if myDate in dataMatrixDict:
            actualDepths = sorted(dataMatrixDict[myDate].keys())

            for index in range(len(actualDepths)-1):
                loDepth = actualDepths[index]
                loTemp = dataMatrixDict[myDate][loDepth][0]
                hiDepth = actualDepths[index+1]
                hiTemp = dataMatrixDict[myDate][hiDepth][0]

                m = findSlope(hiDepth,loDepth,hiTemp,loTemp)
                allDepths = np.arange(loDepth, hiDepth, resolution) #this could be the resolution parameter, we could have a default

                for depth in allDepths:
                    if (depth != loDepth) and (depth != hiDepth):
                        interTemp = findPoint(m,hiDepth,hiTemp,depth)
                        dataMatrixDict[myDate][round(depth,3)] = [interTemp, '']

        myDate += datetime.timedelta(hours=1) #TODO: replace with actualDates, all the dates between my end and start date
        
        #TODO: add verbose parameter that says how many you skipped
        #TODO: write to file, so next time they don't have to run this again
        #TODO: should we give the interpolated data points a flag?
    if verbose:
        print findGaps(dataMatrixDict, startDate, endDate) #TODO: when you have CLI running, you should make a helper function for printing errors/messages
    return dataMatrixDict     

def sunriseToDateTime(sunriseData):
    """Converts sunrise/sunset data into datetime format.
        Args:
            sunriseData: matrix of sunrise/sunset data for year in string format
        Returns:
            matrix of sunrise/sunset data for year in datetime format
    """
    dates = []
    year = sunriseData[:1][0]
    for row in sunriseData[1:]:
        for col in range(len(row)):
            if col == 0:
                day = row[col]
            else:
                month = col 
                if len(row[col]) != 0:
                    pair = []
                    for i in row[col]:
                        if len(str(i)) == 3:
                            i = ''.join(('0',str(i)))
                        pair.append(str(i))
                    start = pair[0]
                    sHour = int(start[:2])
                    sMin = int(start[2:])

                    end = pair[1]
                    eHour = int(end[:2])
                    eMin = int(end[2:])
                    
                    dates.append([datetime.datetime(year, month, day, sHour, sMin),
                                  datetime.datetime(year, month, day, eHour, eMin)])
    return sorted(dates)

## helper functions ##

def interpolate(hiDepth,loDepth,hiTemp,loTemp):
    """Finds interpolated temperature between two known points.

    Args: 
        hiDepth: the known depth higher
        loDepth: the known depth lower
        hiTemp: the known temp at the hiDepth
        loTemp: the known temp at the loDepth

    Returns:
        interpolated temperature rounded to 3 decimal places
    """
    m = findSlope(hiDepth,loDepth,hiTemp,loTemp)
    interpTemp = findPoint(m,hiDepth,hiTemp,targetDepth)

    return round(interpTemp,3)

def findSlope(x1,x2,y1,y2):
    """ Finds the slope of a line.

    Args:
        x1, y1: the coordinates of the first point
        x2, y2: the coordinates of the second point

    Returns: 
        m: the slope of the line 
    """
    m = (y1-y2)/(x1-x2)

    return m

def findPoint(m, x, y, x1):
    """ Given the slope and one point, finds the y value for a second point on the line.

    Args:
        m: the slope of the line
        x, y: the coordinates of a point on the line
        x1: the x value of the point of interest

    Returns: 
        y1: the y value of the point of interest 
    """
    y1 = y - (m * (x-x1))
    y1 = round(y1,3)

    return y1

def formatDateHourly(dateString):
    """ Converts hourly string (in form %Y-%m-%d %H:%M:%S) to a datetime object. 

    Args:   
        dateString: hourly date string

    Returns:
        datetime object
    """
    formattedDate = datetime.datetime.strptime(dateString, '%Y-%m-%d %H:%M:%S')
    
    return formattedDate

def formateDateDaily(dateString):
    """ Converts string (in form %Y-%m-%d) to a datetime object.

    Args:
        dateString: date string
    Returns: datetime object 

    """
    formattedDate = datetime.datetime.strptime(dateString, '%Y-%m-%d ')

    return formattedDate

def dateString(date, year=False):
    """ Converts year/date into a start and end datetime objects. 

    Args:
        date: date string
        year: boolean representing whether the date is a single year or a year/month/day string.
    
    Returns:
        start, end: the start and end datetime objects
    """
    dateString = str(date)
    if year:
        start = formatDateHourly(dateString+'-01-01 00:00:00') #this can be modified to isolate a season
        end = formatDateHourly(dateString+'-12-31 23:00:00')
    else:
        start = formatDateHourly(dateString+' 00:00:00')
        end = formatDateHourly(dateString+' 23:00:00')

    return start, end

def getDates(startYear, endYear):
    """Given a start year and end year, returns a list containing all years in between.

    Args:
        startYear: the starting year
        endYear: the ending year

    Returns:
        dates: list of all years inclusive
    """
    if startYear != endYear:
        dates = range(startYear,endYear)
    else:
        dates = [startYear]

    return dates

def writeCSV(dictionary, nameString):
    """Writes a dictionary into a csv (excel) file.
    
    Args:
        dictionary: the dictionary to be written
        nameString: the name of the file created
    
    Returns:
        writes a file to <nameString>.csv in the current directory.
    """
    w = csv.writer(open(nameString+".csv","w"))
    for date in dictionary.keys():
        for depth in dictionary[date]:
            ls = dictionary[date][depth]
            if isinstance(ls, (list,)):   
                temp = dictionary[date][depth][0]
                flag = dictionary[date][depth][1]
                if len(ls) == 3:
                    frequency = dictionary[date][depth][2]
                    w.writerow([date, depth, temp, flag, frequency])
                else:
                    w.writerow([date, depth, temp, flag])
            else:
                w.writerow([date, depth, dictionary[date][depth]])
