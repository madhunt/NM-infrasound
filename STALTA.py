#!/usr/bin/env python3
"""
Pick events using STALTA (short term average long term average)
Created on 31 Jan 2020

@author: madhunt
"""

import numpy as np
import matplotlib.pyplot as plt
import obspy as obs
import obspy.signal.trigger as tg
import obspy.signal.filter as fl
from os import listdir
from os.path import isfile, join


def procData(data, dec, detrend, filt, **filt_kwargs):
    '''
    Decimates, detrends, and filters data.
    INPUTS:
        data: seismic or infrasound data; type: trace
        dec: factor to downsample the data by; type: integer
        detrend: 'simple','linear','constant','polynomial','spline'; type: str
        filt: type of filter; 'bandpass','lowpass','highpass'; type: str
        filt_kwargs: options depending on filter chosen (eg. freqmin,freqmax); NOTE: df not needed as input
    RETURNS: 
        data_proc: array of processed data
    '''
    # decimate data (decrease sampling rate)
    data_proc = data.copy()
    data_proc.decimate(dec)
    
    # center data at 0
    data_proc.detrend(detrend)
    
    # filter data
    data_proc.filter(filt,**filt_kwargs)
    
    return data_proc


def stalta(data, w_s):
    '''
    Computes STA/LTA ratios from data.
    INPUTS:
        data: seismic or infrasound data (preferably processed); type: trace
        w_s: [STA,LTA] window sizes in seconds; type: 2 value array
    RETURNS: 
        cft: array of STALTA ratios
        samprate: sample rate (if using processed data, this is orig_samprate/dec)
    '''
    
    # STA and LTA window lengths in terms of data points
    samprate = data.stats.sampling_rate
    STA_w = (samprate * w_s[0])
    LTA_w = (samprate * w_s[1])
    
    # calculate STA/LTA ratios
    cft = tg.classic_sta_lta(data,STA_w,LTA_w)
    
    return cft,samprate


def pickEvents(path_data, w_s, thres, dec, detrend, filt, **filt_kwargs):
    '''
    Loads in data, runs stalta function, and returns list of events
    INPUTS:
        path_data: path to data (stored in mseed format); type: str
        w_s: [STA,LTA] window sizes in seconds; type: 2 entry array
        thres: [lower, upper] threshold where trigger is deactivated or activated; type: 2 entry array
        dec: factor to downsample the data by; type: int
        detrend: 'simple','linear','constant','polynomial','spline'; type: str
        filt: type of filter; 'bandpass','lowpass','highpass', etc.; type: str
        filt_kwargs: options depending on filter chosen; type: str
            if filt='bandpass': freqmin=lower cutoff frequency, freqmax=upper cutoff
            if filt='lowpass': freq=cutoff frequency
            if filt='highpass': freq=cutoff frequency
    RETURNS:
        events: list with [date/filename (str), event times (int)]
    '''
    # make a list of data files in given directory
    files = [f for f in listdir(path_data) if isfile(join(path_data, f))]
    
    # make empty array to store [filename, event time]
    events = [[],[]] 
    
    # loop through each file
    for filenum in range(0,np.size(files)):
        print('File Number '+str(filenum)+' Out of '+str(np.size(files)))
        
        # read in the data
        mseed = obs.read(join(path_data, files[filenum]))
        data = mseed[0]
        
        # process data
        data_proc = procData(data, dec, detrend, filt, **filt_kwargs)
        
        # calculate STALTA ratios
        cft,samprate = stalta(data_proc, w_s)
        
        npoints = 3600 * samprate
        
        ############################################
        #TODO figure out why some files are too big or too small
        #TODO find a way to fix this and get events from these files
        #HACK for now, just skip these files
        if np.size(cft) < npoints:
            print('uh oh... file number ' + str(filenum) + 'is too small!')
            continue
#        # if there is not a full hour of data, pad end with zeros
#        if np.size(cft) < npoints:
#            npad = npoints-np.size(cft)
#            cft = np.pad(cft, npad, 'constant')
        
        if np.size(cft) > npoints:
            print('uh oh... file number ' + str(filenum) + 'is too big???')
            continue
        ############################################ 
        
        # choose events
        new_events = tg.trigger_onset(cft,thres[1],thres[0])
        
        # if there are new events, store them
        if len(new_events) != 0:
            # store each event in the same hour in a separate row
            for n in range(len(new_events)):
                # store hour / filename
                events[0].append(files[filenum]) 
                
                # store the event time 
                events[1].append(int(new_events[n][0]))
    
    return events






def runFuncs():
    # path to data (mseed)
    path_data = '/home/mad/Documents/Research2020/MSEEDdata/WCYN/'
    # no final / so MMTN is in filename
    path_save = '/home/mad/Documents/Research2020/picks/WCYN' 
    
    # [STA window length (sec), LTA window length (sec)]
    w_s = [5,50]
    # [lower, upper threshold values where trigger is deactivated, activated]
    thres = [0.5,7]
    
    dec = 5
    detrend = 'linear'
    filt = 'bandpass'
    
    # for filtering with bandpass
    freqmin = 1
    freqmax = 7
    
    events = pickEvents(path_data, w_s, thres, dec, detrend, filt, freqmin=freqmin,freqmax=freqmax)
    
    # save events to computer
    np.save((path_save + 'events'), events)

    
    return








