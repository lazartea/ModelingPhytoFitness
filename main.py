import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import random
import csv

from movementStrategies import *
from createFigures import *
from statistics import *
from formatData import *

def main(filename, singleYear, hourly):
    dataMatrix = createDataMatrix(filename,'hourly')    
    dictionary = createDictionary(dataMatrix)
    
    start, end = dateString(singleYear, hourly)
    
    dateDicts = groupAllDates(dictionary)
    averageDict = averageOverYears(dateDicts)
    dataMatrixDict = fillGapsInData(dictionary, averageDict, singleYear)
    
    dataDict = extendDataMatrix(dataMatrixDict, start, end)
    sunriseData = createSunriseData('sunrise/'+str(singleYear)+'sunrise.txt', singleYear)
    sunriseData = sunriseToDateTime(sunriseData)
    fitnessDict = createFitnessDict(dataDict, sunriseData, start, end)
    
    return dataDict, fitnessDict, sunsetMatrix

data, fitness, sunset = main('sensorsparklinglakewatertemphourly.txt', 2005, True) #use these to call movement functions
#ex: 
#start, end = dateString(2000, hourly) #where 2000 is year of interest
#circadianMovement(data, fitness, sunset, start, end, 'slow')
