#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 2020

@author: madhunt
"""
import numpy as np

def geoMid(lat,long):
    '''
    Find geographic midpoint.
    Using method described at http://www.geomidpoint.com/calculation.html
    INPUTS:
        lat: array of latitudes in degrees
        long: array of longitudes in degrees
    RETURNS:
        latAvg, longAvg: Geographic midpoint in degrees; type: float
    '''
    # convert to radians
    latr,longr = map(np.radians, [lat,long])
    
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


def avglatlong(path):
    '''
    Calculates average latitude and longitude using geometric midpoint.
    INPUTS:
        path: path to position file; type: str
    RETURNS:
        latAvg,longAvg,zAvg: average latitude, longitude, and elevation; type: float
    '''
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
    
    latAvg,longAvg = geoMid(lat,long)
    zAvg = np.average(z)
    
    return latAvg,longAvg,zAvg


def geoDist(lat1,long1,z1,lat2,long2,z2):
    '''
    Calculates distance in meters between two points on the Earth using the Haversine formula.
    Using method described on https://www.movable-type.co.uk/scripts/latlong.html
    INPUTS:
        lat1,long1,z1: latitude, longitude, and elevation of point 1 (lat/long in decimal notation)
        lat2,long2,z2: latitude, longitude, and elevation of point 2 (lat/long in decimal notation)
    RETURNS:
        dist: straight-line distance between the points; type: float
    '''
    # convert to radians
    long1, lat1, long2, lat2 = map(np.radians, [long1, lat1, long2, lat2])
    
    # define radius of earth in m
    R = 6371 * 10**3
    
    # calculate some differences
    dlat = lat2 - lat1
    dlong = long2 - long1
    dz = z2 - z1
    
    # use Haversine formula
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlong/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    xydist = R * c
    
    # compute distance including elevation
    dist = np.sqrt(xydist**2 + dz**2)
    
    return dist


def approxTimeBtw(pathM,pathW,vs):
    '''
    INPUTS:
        pathM: path to gps position file for station M; type: str
        pathW:
        vs: speed of sound (m/s), type: float
    RETURNS:
        tapprox: approximate time (s) for sound to travel between stations; type: float
    '''
    # calculate average lat/long of two stations
    latM,longM,zM = avglatlong(pathM)
    latW,longW,zW = avglatlong(pathW)
    
    # calculate distance between stations
    dist = geoDist(latM,longM,zM,latW,longW,zW)
            
    # approximate travel time in sec between stations
    t_btw = dist/vs
    
    return t_btw


def main():
    # paths to position data
    pathM = '/home/mad/Documents/Research2020/gps/MMTNlog/pos.sta1'
    pathW = '/home/mad/Documents/Research2020/gps/WCYNlog/pos.sta2'
    
    # use approx speed of sound in m/s
    vs = 343
    
    # predicted time between stations
    tpred = approxTimeBtw(pathM,pathW,vs)
    
    return tpred


#main()
