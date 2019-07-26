import math 
import numpy as np
import decimal
import datetime
import random
import sys 
## movement patterns ##

def fitnessFunction(temp, date):
    #these parameters come from Colin Kramer, for Cyanobacteria Synechococcus (which are in high abundance in Sparkling Lake)
    b1 = 14.59829888 #birth rate at 0C
    b2 = 0.008383057 #change in birth rate due to T
    d0 = 9.301242586 #temperature-independent mortality constant
    d1 = 5.545155413 #exp. changes to mortality due to incr. T
    d2 = 0.016661372 #exp. changes to mortality due to incr. T
    tOpt = 33.95619378 #the temperature at which growth rate is maximized
    growth = (b1 * math.exp(b2*temp)) - (d0 + ((b1*b2)/d2) * math.exp((b2-d2)*tOpt) * math.exp(d2 * temp))
    
    return growth

def doubleExponentialGrowthRate(temp):
    #original version of fitness function used in thesis
    b1 = 1.174 #birth rate at 0C
    b2 = 0.064 #change in birth rate due to T
    d0 = 1.119 #temperature-independent mortality constant
    d1 = 0.267 #exp. changes to mortality due to incr. T
    d2 = 0.103 #exp. changes to mortality due to incr. T

    growth = (b1 * math.exp(b2 * temp) - (d0 + d1 * math.exp(d2 * temp)))/1.2
    return growth 

def createFitnessDict(dataMatrixDict, sunriseDict, start, end):
    """Creates dictionary of fitnesses across entire lake for quicker processing.

    Args: 
        dataMatrixDict: the dictionary of depths and temperatures
        start: the start date
        end: the end date
    Returns:
        fitnessDict: dictionary of fitnesses at each point in time and space
    """
    fitnessDict = {}
    for date in sorted(dataMatrixDict.keys()):
        #sunrise, sunset = findSunsetTimes(date, sunriseDict)
        if date >= start and date <= end:
            if date not in fitnessDict: 
                fitnessDict[date] = {}
            for depth in dataMatrixDict[date].keys():
                if depth not in fitnessDict[date]:
                    fitnessDict[date][depth] = fitnessFunction(dataMatrixDict[date][depth][0], date)
    
    return fitnessDict

def fitnessAtStableDepth(dataMatrixDict, fitnessDict, targetDepth, start, end, hourly):
    """Finds fitness at each point in time at a specific depth.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        targetDepth: the depth of interest
        start: the start date
        end: the end date
    Returns: 
        dateList = list of each date where fitness was calculated
        fitnessList = each fitness found 
    """
    dateList = sorted(getCommonElements(getDatesBetweenRange(start, end, hourly), dataMatrixDict.keys()))
    tempList, fitnessList = [], []
    for date in dateList:
        for depth in dataMatrixDict[date].keys():
            if depth == targetDepth:
                temp = dataMatrixDict[date][depth][0]
                fitnessList.append(fitnessDict[date][depth])

    return dateList, tempList, fitnessList 

def circadianMovement(dataMatrixDict, fitnessDict, sunsetMatrix, start, end, speed, seed=False):
    """Simulates Diel Vertical Migration, finds fitness.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        sunsetMatrix: the matrix of sunrise and set times for each day
        start: the start date
        end: the end date
        speed: the speed of movement, representing upper and lower bounds of observed speeds
    Returns: 
        dateList, depthList, tempList, fitnessList
    """
    if speed == 'fast':
        distance = 4
    elif speed == 'slow':
        distance = 2
    else:
        print "Error: incorrect parameters" 
        sys.exit()
    if seed: 
        random.seed(2)

    dateList = sorted(getCommonElements(getDatesBetweenRange(start, end, True), dataMatrixDict.keys()))
    depthList, tempList, fitnessList = [], [], []
    
    for date in dateList:
        sunrise, sunset = findSunsetTimes(date, sunsetMatrix)
        if len(depthList) == 0:
            depth = random.choice(dataMatrixDict[date].keys())
            previousLocation = depth
            temp = dataMatrixDict[date][previousLocation][0]

        else:
            depthsNum = range(len(dataMatrixDict[date].keys()))
            allDepths = sorted(dataMatrixDict[date].keys())
            if previousLocation in allDepths:
                index = allDepths.index(previousLocation) 
            else:
                #TODO: break (count how often this happens)
                #TODO: separate checking for jumps and handling hitting bottom/surface
                previousLocation = min(allDepths, key=lambda x:abs(x-previousLocation))
                index = allDepths.index(previousLocation) 

            if date >= sunrise and date < sunset:
                index = index - distance
            else:
                index = index + distance
            
            if index in depthsNum:
                previousLocation = allDepths[index]
        depthList.append(previousLocation)
        tempList.append(dataMatrixDict[date][previousLocation][0])
        fitnessList.append(fitnessDict[date][previousLocation])

    return dateList, depthList, tempList, fitnessList

def hillClimbingMovement(dataMatrixDict, fitnessDict, start, end, hourly, seed=False):
    """Hill Climbing movement pattern.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        start: the start date
        end: the end date
        seed: optional. used to set a seed for testing/repeatability purposes
        numVers: optional. used to diffrentiate when running multiple times
    Returns: 
        dateList, depthList, tempList, fitnessList

    """
    dateList = sorted(getCommonElements(getDatesBetweenRange(start, end, hourly), dataMatrixDict.keys()))

    depthList, tempList, fitnessList = [], [], []

    for date in dateList:
        if len(depthList) == 0:
            if seed:
                random.seed(2)
            depth = random.choice(dataMatrixDict[date].keys())
            temp = dataMatrixDict[date][depth][0]
            fitness = fitnessDict[date][depth]
            previousLocation = depth
                                  
        else: 
            depthsNum = range(len(dataMatrixDict[date].keys()))
            allDepths = sorted(dataMatrixDict[date].keys())
            if previousLocation in allDepths:
                index = allDepths.index(previousLocation) 
            else:
                previousLocation = min(allDepths, key=lambda x:abs(x-previousLocation))
                index = allDepths.index(previousLocation)
            
            comparisons = []
            curTemp = dataMatrixDict[date].get(previousLocation, False)
            if curTemp:
                temp = curTemp[0] 
                fitness = fitnessDict[date].get(previousLocation, False)

            if (index+1) in depthsNum:
                hiDepth = allDepths[index + 1]
                hiTemp = dataMatrixDict[date].get(hiDepth, False)

                if hiTemp:
                    hiTemp = hiTemp[0]
                    hiFitness = fitnessDict[date].get(hiDepth, False)
                    if hiFitness > fitness:
                        previousLocation = hiDepth
                        temp = hiTemp
                        fitness = hiFitness

            if (index-1) in depthsNum:
                loDepth = allDepths[index - 1]
                loTemp = dataMatrixDict[date].get(loDepth, False)

                if loTemp:
                    loTemp = loTemp[0]
                    loFitness = fitnessDict[date].get(loDepth, False)
                    if loFitness > fitness:
                        previousLocation = loDepth
                        temp = loTemp
                        fitness = loFitness
              
        depthList.append(previousLocation)
        tempList.append(temp)
        fitnessList.append(fitness) 

    return dateList, depthList, tempList, fitnessList

def oracleMovement(dataMatrixDict, fitnessDict, start, end, hourly): 
    """The Oracle finds the best depth (by maximizing fitness) at every point in time.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        start: the start date
        end: the end date
    Returns: 
        dateList = list of each date where fitness was calculated
        fitnessList = each fitness found 
    """
    dateList = sorted(getCommonElements(getDatesBetweenRange(start, end, hourly), dataMatrixDict.keys()))
    depthList, fitnessList, tempList = [], [], []

    for date in dateList:
        previousFitness = -100
        #fitness = [doubleExponentialGrowthRate(d[0]) for depth in dataMatrixDict[date].keys()] list comprehension
        for depth in dataMatrixDict[date].keys():
            fitness = fitnessDict[date][depth]
            if fitness > previousFitness:
                previousFitness = fitness
                previousLocation = depth 
                temp = dataMatrixDict[date][depth][0]

        depthList.append(previousLocation)
        fitnessList.append(previousFitness)
        tempList.append(temp)

    return dateList, depthList, tempList, fitnessList

def randomWalk(dataMatrixDict, fitnessDict, start, end, hourly):
    """Moves to a random location at each point in time, calculates fitness.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        start: the start date
        end: the end date
        hourly: boolean representing the resolution of the data (hourly or daily)
    Returns: 
       dateList, depthList, tempList, fitnessList
    """
    dateList = getCommonElements(getDatesBetweenRange(start, end, hourly), dataMatrixDict.keys())
    depthList, fitnessList, tempList = [], [], []
    
    for date in dateList:
        depth = random.choice(dataMatrixDict[date].keys())
        temp = dataMatrixDict[date][depth][0]
        fitness = fitnessDict[date][depth]
        depthList.append(depth)
        fitnessList.append(fitness)
        tempList.append(temp)

    return dateList, depthList, tempList, fitnessList

def randomWalkDirectional(dataMatrixDict, fitnessDict, start, end, probabilityFactor, hourly):
    """Moves one step up or down (using randomness) at each interval, calculates fitness.

    Args:
        dataMatrixDict: the dictionary of depths and temperatures
        fitnessDict: dictionary of fitnesses at each point in time and space
        start: the start date
        end: the end date
        probabilityFactor: the weight of the coin flipped
        hourly: boolean representing the resolution of the data (hourly or daily)
    Returns: 
       dateList, depthList, tempList, fitnessList
    """
    dateList = getCommonElements(getDatesBetweenRange(start, end, hourly), dataMatrixDict.keys())
    depthList, tempList, fitnessList = [],[],[]
    
    for date in dateList:
        depths = sorted(dataMatrixDict[date].keys())
        indices = range(len(depths)-1)
        if len(depthList) == 0:
            depth = random.choice(depths)
            previousIndex = depths.index(depth)
            
        else:
            rand = random.random()
            if rand > probabilityFactor:
                if previousIndex - 1 >= 0 and previousIndex-1 in indices:
                    previousIndex = previousIndex - 1
                else:
                    if previousIndex not in indices:
                        previousIndex = min(indices, key=lambda x:abs(x-previousIndex))
                #else, previousIndex stays the same because we can't decrease from 0 (at the surface)
            else:
                if previousIndex + 1 < len(depths):
                    previousIndex = previousIndex + 1
                else: 
                    previousIndex = len(depths) - 1 #for most cases, previousIndex will stay the same. but if the max lake depth has changed, reset to deepest 
        depth = depths[previousIndex]
        temp = dataMatrixDict[date][depth][0]
        fitness = fitnessDict[date][depth]
        previousIndex = depths.index(depth) 
        depthList.append(depth)
        fitnessList.append(fitness)
        tempList.append(temp)

    return dateList, depthList, tempList, fitnessList


## helper functions ##
def getCommonElements(listA, listB):
    """ Takes in two lists and returns a list with their common elements.

    Args:
        listA, listB: two lists

    Returns:
        common elements of the lists
    """
    setA = set(listA)
    setB = set(listB)
    
    return list(setA.intersection(setB))

def findSunsetTimes(currentDate, sunsetMatrix):
    """Finds sunrise and sunset times for the given day.

    Args:
        currentDate: the date of interest
        sunsetMatrix: the matrix of all sunrise/sunset data

    Returns:
        sunrise time, sunset time 

    """
    for row in sunsetMatrix:
        if row[0].date() == currentDate.date():
            return row[0], row[1]

def getDatesBetweenRange(start, end, hourly):
    """Returns all dates between start and end date.
    Args:
        start: the start date
        end: the end date
        hourly: boolean representing the resolution of the data (hourly or daily)
    Returns: 
       dateList
    """
    dateList = []
    if hourly:
        currentDate = start
        dateList.append(start)
        while currentDate < end:
            currentDate = currentDate + datetime.timedelta(hours = 1)
            dateList.append(currentDate)

    else:
        currentDate = start
        dateList.append(start)
        while currentDate < end:
            currentDate = currentDate + datetime.timedelta(days = 1)
            dateList.append(currentDate)

    return dateList