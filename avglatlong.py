#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get latitude and longitude from station logs and find averages
Created on Mon Feb 17 2020

@author: madhunt
"""

import numpy as np

def avglatlong(path):
    '''
    Get latitude, longitude, and elevation; convert to decimal notation; calculate the geographic midpoint (average)
    INPUTS:
        path: path to position file; type: str
    RETURNS:
        latAvg,longAvg,zAvg: average latitude, longitude, and elevation; type: float
    '''
    
    # open file
    file = open(path,mode='r')
    
    # make list with position strings
    lines = file.readlines()
    
    # get lat, long, and elevation
    pos_str = [row.split(' ')[3:6] for row in lines]
    
    def degtodec(di,deg,m,s):
        '''
        Convert degrees/minutes/seconds to decimal
        INPUTS:
            di: compass direction ('N','S','E','W'); type: str
            deg: degrees; type: int
            m: minutes; type: int
            s: seconds; type: int or float
        RETURNS:
            l: lat or long in degrees; type: float
        '''
        l = float(deg) + float(m)/60 + float(s)/3600
        
        if di == 'W' or di == 'S':
            l *= -1
        
        return l
    
    # convert lat long to decimal
    pos = [[],[],[]]
    
    for row in pos_str:
        # make elevation into an integer
        pos[2].append(int(row[2][1:-2]))
        
        for n in range(2):
            # convert lat and long into degrees
            deg,m,s = row[n].split(':')
            di = deg[0]
            deg = deg[1:]
            l = degtodec(di,deg,m,s)
            
            pos[n].append(l)
    
    # convert to numpy array because now it's math time
    pos = np.asarray(pos)
    
    lat = pos[0]
    long = pos[1]
    z = pos[2]
    
    def geoMid(lat,long):
        '''
        Find the Geographic midpoint (average assuming perfect sphere)
        Using method described on http://www.geomidpoint.com/calculation.html
        INPUTS:
            lat: array of latitudes in degrees
            long: array of longitudes in degrees
        RETURNS:
            latAvg, longAvg: Geographic midpoint in degrees; type: float
        '''
        # convert to radians
        latr = lat * np.pi/180
        longr = long * np.pi/180
        
        # convert to cartesian points
        x = np.cos(latr) * np.cos(longr)
        y = np.cos(latr) * np.sin(longr)
        z = np.sin(latr)
        
        # compute average
        xAvg = np.average(x)
        yAvg = np.average(y)
        zAvg = np.average(z)
        
        # convert back to lat/long
        latAvgr = np.arctan2(zAvg,np.sqrt(xAvg**2 + yAvg**2))
        longAvgr = np.arctan2(yAvg,xAvg)
        
        # convert back to degrees
        latAvg = latAvgr * 180/np.pi
        longAvg = longAvgr * 180/np.pi
        
        return latAvg,longAvg
    
    latAvg,longAvg = geoMid(lat,long)
    zAvg = np.average(z)
    
    return latAvg,longAvg,zAvg


# get position files
pathM = '/home/mad/Documents/Research2020/gps/MMTNlog/pos.sta1'
pathW = '/home/mad/Documents/Research2020/gps/WCYNlog/pos.sta2'


latM,longM,zM = avglatlong(pathM)
latW,longW,zW = avglatlong(pathW)

print(latM,longM)
print(latW,longW)

# can then put these coordinates in Google Earth


