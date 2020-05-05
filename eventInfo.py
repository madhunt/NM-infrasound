#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 2020

@author: madhunt
"""

import numpy as np
import matplotlib.pyplot as plt

import csv


def readData(pathTrig):
    
    trig = np.load(pathTrig,allow_pickle='TRUE').item()
    
    return trig


def sigToNoise(data):
    '''
    Computes signal-to-noise ratio of data, assuming data is centered at 0.
    INPUTS:
        data; data with signal triggered at len/2; type: list or array
    RETURNS:
        s2n: signal-to-noise ratio; type: float
    '''
    n = len(data)
    
    Snoise = np.var(data[0:int(0.5*n)])
    Sboth = np.var(data[int(0.5*n):])
    Ssignal = Sboth - Snoise
    
    s2n = Ssignal/Snoise
    
    return s2n


def RMS(data):
    '''
    Computes RMS amplitude.
    '''
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
            s1 = sigToNoise(data1)
            s2 = sigToNoise(data2)
            database['sigToNoise'].append((s1,s2))
            
            # find variance
            LTA1 = 30/120/2
            LTA2 = 1-LTA1
            
            v1 = np.var(data1[int(LTA1*n1) : int(LTA2*n1)])
            v2 = np.var(data2[int(LTA1*n2) : int(LTA2*n2)])
            database['variance'].append((v1,v2))
            
            # find RMS amplitude
            STA1 = 2/120/2
            STA2 = 1-STA1
            
            rms1 = RMS(data1[int(STA1*n1) : int(STA2*n1)])
            rms2 = RMS(data2[int(STA1*n2) : int(STA2*n2)])
            database['amp_RMS'].append((rms1,rms2))
            
            # find peak to peak amplitude          
            
            # increase event ID
            ID += 1
    
    return database
    

def plot1event(ID,trig):
    # loop through days and events, etc. to find event matching given ID
    for day in trig:
        events = trig.get(day)
        for event in events:
            
            if ID == event['ID']:
                data1 = event['data1']
                data2 = event['data2']
                
                n1 = len(data1)
                n2 = len(data2)
                
                time = event.get('time')
                timeDiff = time[1] - time[0]
                
                timeaxis1 = np.linspace(-60,60,int(n1/2))
                timeaxis2 = timeaxis1 + timeDiff

                plt.close('all')
                plt.figure()
                
                plt.subplot(211)
                plt.plot(timeaxis1,data1[int(0.25*n1):int(0.75*n1)])
                plt.title(event['stations'][0] + ': ' + time[0].strftime("%d %B %Y, %H:%M:%S"))
                plt.xlim([-30,30])
                plt.ylabel('Amplitude (Pa)')
                plt.text(0, 0,'1 to 7Hz Filter',transform=plt.gca().transAxes)
                
                plt.subplot(212)
                plt.plot(timeaxis2,data2[int(0.25*n2):int(0.75*n2)])
                plt.title(event['stations'][1] + ': '  + time[1].strftime("%d %B %Y, %H:%M:%S"))
                plt.xlim([-30,30])
                plt.ylabel('Amplitude (Pa)')
                plt.xlabel('Time from Event (sec)')
                plt.text(0, 0,'1 to 7Hz Filter',transform=plt.gca().transAxes)
                
                plt.suptitle('Event ' + str(ID))
    return


def saveDatabase(my_dict, path):
    with open(path, "wb") as f_output:
        csv_output = csv.writer(f_output)
        csv_output.writerow(['ID', 'timeDiff', 'sigToNoise', 'variance', 'amp_RMS', 'amp_p2p'])
    
        for key in sorted(my_dict.keys()):
            csv_output.writerow([key] + my_dict[key])
            
    return


def main():
    # path to event triggers and to save data
    path = '/home/mad/Documents/Research2020/trig/'
    
    trig = readData(path + 'eventTrig_subset.npy')
    database = properties(trig)
    
    # event ID to plot
    ID = 42
        
    plot1event(ID,trig)
    
    #saveDatabase(database,(path+'database_subset.csv'))
    
    return trig, database


trig,database = main()
