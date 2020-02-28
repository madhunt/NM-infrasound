#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to plot events.

Created on Wed Feb 26 2020

@author: madhunt
"""

import numpy as np
import obspy as obs
from os import listdir
from os.path import join
import matplotlib.pyplot as plt
from STALTA import procData # import procData function from STALTA.py, by madhunt


def plotEvent(path,events,eventNum,plot_window,**proc_kwargs):
    '''
    INPUTS:
        path: path to raw mseed data; type: str
        events: list of events in format [date/filename (str), event times (int)]
        eventNum: index of one event in events list to plot
        plot_window: [time before event, time after event] to include in image (in s)
        proc_kwargs: options to pass to procData()
            (dec, detrend, filt, freqmin, freqmax)
    RETURNS:
        plots the event in new figure
    '''
    # get event
    datestr = events[0][eventNum]
    n_arr = int(events[1][eventNum])
    
    # open unprocessed data file
    mseed = obs.read(join(path,datestr))
    data = mseed[0]
    
    # process data (to plot)
    data_proc = procData(data=data,**proc_kwargs)
    samprate = data_proc.stats.sampling_rate
    
    # find indices to plot between
    n_start = int(n_arr - plot_window[0] * samprate)
    n_stop = int(n_arr + plot_window[1] * samprate)
    
    # if final index is too large, load in next consecutive file
    n_tot = np.size(data_proc)
    if n_stop > n_tot:
        # load in next file
        i = listdir(path).index(datestr)
        
        assert i+1 >= len(listdir(path)), ('Event occurs in last file of data')
        
        datestr2 = listdir(path)[i+1]
        
        data2 = obs.read(join(path,datestr2))[0]
        data2_proc = procData(data=data2,**proc_kwargs)
        
        # add total indices in data to data2
        data2_proc = [x + n_tot for x in data2_proc]
        
        # now concatinate together so current data is long enough
        data_proc = data_proc + data2_proc
    
    # now we have enough data to plot
    data_plot = data_proc[n_start:n_stop]
    
    # make an array of times
    t_arr = n_arr / samprate
    t_start = t_arr - plot_window[0]
    t_stop = t_arr + plot_window[1]
    t_plot = np.linspace(t_start,t_stop,np.size(data_plot))
    
    # plot the figure
    fig,ax = plt.subplots()
    plt.plot(t_plot,data_plot,'k')
    plt.xlabel('Time (s)')
    ##TODO: units for y-axis (Pa?)
    plt.ylabel('Amplitude (Pa?)')
    plt.title(datestr)
    # add text to indicate that data was filtered
    plt.text(0,0,'Bandpass Filter: '+str(freqmin)+'Hz to '+str(freqmax)+'Hz',ha='left',va='bottom',transform=ax.transAxes)
    plt.show()
    
    return




# load in events
eventsW = np.load('/home/mad/Documents/Research2020/picks/WCYNevents.npy')
# path to original mseed data
path = '/home/mad/Documents/Research2020/MSEEDdata/WCYN/'

dec = 5
detrend = 'linear'
filt = 'bandpass'
freqmin = 1
freqmax = 7
eventNum = 14
plot_window = [30,30]

plotEvent(path,eventsW,eventNum,plot_window, dec=dec, detrend=detrend, filt=filt, freqmin=freqmin,freqmax=freqmax)





