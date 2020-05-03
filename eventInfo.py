#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 2020

@author: madhunt
"""

import numpy as np
import matplotlib.pyplot as plt

import datetime

import csv


def readData(pathTrig):
    trig = np.load(pathTrig,allow_pickle='TRUE').item()
    return trig


def sigToNoise(data):
    '''
    Computes signal-to-noise ratio of data; from scipy.stats v0.14.0
    INPUTS:
        data; type: list or array
    RETURNS:
        s2n: signal-to-noise ratio; type: float
    '''
    data = np.asanyarray(data)
    mean = data.mean()
    std = data.std()
    
    s2n = np.where(std == 0, 0, mean/std)
    
    return float(s2n)


def RMS(data):
    data = np.asarray(data)
    
    rms = np.sqrt(np.mean(data**2))
      
    return rms


def properties(trig):
    '''
    Creates dictionary with various properties of events given in trig.
    INPUTS:
        trig: dictionary with events, found in STALTA.py
    RETURNS:
        database: dictionary with various properties
            'ID': event number; type: int
            'timeDiff': difference between arrival times at 2 stations, in sec
            'sigToNoise': signal-to-noise ratio; type: tuple
            'variance': variance of signal; type: tuple
            'amp_RMS': RMS amplitude of signal at each station; type: tuple
            'amp_p2p': peak-to-peak amplitude; type: tuple
            
    '''
    database = {}   
    database['ID'] = []
    database['timeDiff'] = []
    database['sigToNoise'] = []
    database['variance'] = []
    database['amp_RMS'] = []
    database['amp_p2p'] = []    
    
    # event identifier
    ID = 0
    
    # loop through days
    for day in trig:
        
        events = trig.get(day)
        
        # loop through events
        for event in events:
            database['ID'].append(ID)
            
            time = event.get('time')
            data1 = event.get('data1')
            data2 = event.get('data2')
            
            n1 = len(data1)
            n2 = len(data2)
            
            # find difference in arrival times
            timeDiff = time[1] - time[0]
            database['timeDiff'].append(timeDiff)
            
            # find signal to noise ratio 
            s1 = sigToNoise(data1[int(0.25*n1) : int(0.75*n1)])
            s2 = sigToNoise(data2[int(0.25*n2) : int(0.75*n2)])
            database['sigToNoise'].append((s1,s2))
            
            # find variance
            #TODO: should this be the variance across the entire window or just the signal?
            v1 = np.var(data1)
            v2 = np.var(data2)
            database['variance'].append((v1,v2))
            
            # find RMS amplitude
            #TODO: should this be just around the event or the entire window
            rms1 = RMS(data1[int(0.45*n1) : int(0.65*n1)])
            rms2 = RMS(data2[int(0.45*n2) : int(0.65*n2)])
            database['amp_RMS'].append((rms1,rms2))
            
            # find peak to peak amplitude
            #TODO: 
            
            
            
            # increase event ID
            ID += 1
            
    
    return database
    

def plot1(ID,trig):
    
    for day in trig:
        events = trig.get(day)
        for event in events:
            if ID == event['ID']:
                data1 = event['data1']
                data2 = event['data2']
                
                n1 = len(data1)
                n2 = len(data2)
                
                #tevent = event['time']
                #dt = datetime.timedelta(seconds=0,minutes=1)
                
                #time1 = np.linspace(tevent[0]-dt,tevent[0]+dt,n1)
                #time2 = np.linspace(tevent[1]-dt,tevent[1]+dt,n2)
                
                time1 = np.linspace(-1,1,n1)
                time2 = np.linspace(-1,1,n2)
                

                plt.close('all')
                plt.figure()
                
                plt.subplot(211)
                plt.plot(time1,data1)
                
                plt.subplot(212)
                plt.plot(time2,data2)
                
                
                
                
    
    return



def saveDatabase():
    
    return

    
    

def main():
    # path to event triggers
    pathTrig = '/home/mad/Documents/Research2020/trig/eventTrig.npy'
    
    trig = readData(pathTrig)
    database = properties(trig)
    
    
    ID = 14
    
    plot1(ID,trig)
    
    
    return trig, database


trig,database = main()





