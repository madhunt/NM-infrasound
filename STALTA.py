#!/usr/bin/env python3
"""
Created on 31 Jan 2020

@author: madhunt
"""

import numpy as np
import obspy as obs
import obspy.signal.trigger as tg

from obspy.core import Stream

from os import listdir
from os.path import isfile, join

import datetime

import latlong

def procData(data, dec, detrend, filt, **filt_kwargs):
    '''
    Decimates, detrends, and filters data.
    INPUTS:
        data: seismic or infrasound data; type: trace/stream
        dec: factor to downsample the data by; type: integer
        detrend: 'simple','linear','constant','polynomial','spline'; type: str
        filt: type of filter; 'bandpass','lowpass','highpass'; type: str
        filt_kwargs: options depending on filter chosen (eg. freqmin,freqmax); NOTE: df not needed as input
    RETURNS: 
        dataProc: array of processed data
    '''
    # decimate data (decrease sampling rate)
    dataProc = data.copy()
    dataProc.decimate(dec)
    
    # center data at 0
    dataProc.detrend(detrend)
    
    # filter data
    dataProc.filter(filt,**filt_kwargs)
    
    return dataProc


def dayDic(path):
    '''
    Makes a dictionary, grouping by day.
    INPUTS:
        path: path to data; type: str
    RETURNS:
        dic: type: dictionary
    '''
    # make a list of all files in directory
    files = [f for f in listdir(path) if isfile(join(path, f))]
    
    # sort these files chronologically
    files.sort()
    
    # create empty dictionary
    dic = {}
    
    for file in files:
        # find day
        day = int(file.split('.')[1])
        
        # if this is the first instance
        if dic.get(day) is None:
            dic[day] = [file]
        # otherwise, add new hour in the day
        else:
            dic[day].append(file)
    
    return dic


def dataToSave(data, tstart, tstop):
    i=0
    
    for hr in data:
        hrstart = hr.stats['starttime']
        hrstop = hr.stats['endtime']
        
        # if the desired start time is within the hour
        if tstart > hrstart and tstart < hrstop:
            
            # if the stop time is also in that hour
            if tstop < hrstop:    
                istart = int((tstart-hrstart)*hr.stats['sampling_rate'])
                istop = int((tstop-hrstart)*hr.stats['sampling_rate'])
                # save this data
                eventData = hr[istart:istop]
                
            # if stop time is in next hour
            else:                
                istart = int((tstart-hrstart)*hr.stats['sampling_rate'])
                # stop time is in the next hour
                nexthr = data[i+1]
                nexthrstart = nexthr.stats['starttime']
                istop = int((tstop-nexthrstart)*nexthr.stats['sampling_rate'])
                # save this data
                eventData = hr[istart:]
                eventData.append(nexthr[:istop])
        
        else:
            # this is not the hour you are looking for
            continue
        i += 1
        
        return eventData


def pickEvents(path1,path2, w_s, thres, coSumThres, dec, detrend, filt, **filt_kwargs):
    '''
    Loads in data, runs stalta function, and returns list of events using coinidence trigger
    INPUTS:
        path1: path to data for station 1 (stored in mseed format); type: str
        path2: path to data for station 2 (stored in mseed format); type: str
        
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
    
    dic1 = dayDic(path1)
    dic2 = dayDic(path2)
    
    trig = {}
    
    # find predicted time for event to travel between stations
    tpred = latlong.main()
    
    ID = 0
    
    # loop through each day (each key in dictionary 1)
    for day in dic1:
        files1 = dic1.get(day)
        files2 = dic2.get(day)
        
        # if data doesn't exist for one of the stations, ignore the day
        if (files1 is None) or (files2 is None):
            continue
        
        # make empty stream
        data = Stream()
        
        for path in path1,path2:
            if path==path1: files=files1
            if path==path2: files=files2
            
            for file in files:
                # load in the data for one files (one hour)
                st = obs.core.read(join(path,file))
                
                # process (filter, decimate, detrend) the data
                dataProc = procData(st, dec, detrend, filt, **filt_kwargs)
                
                # store processed data
                data += dataProc
            
            if path==path1: data1=data
            if path==path2: data2=data
        
        # use coincidence triggering on processed data from both stations
        trigDay = tg.coincidence_trigger('recstalta',thres[1],thres[0],
                                      data1+data2,coSumThres,sta=w_s[0],lta=w_s[1])
        
        # empty list for REAL events
        realTrigDay = []
        
        # property checking
        for event in trigDay:
            stations = event.get('stations')
            if stations[0] != 'MMTN':
                # 'event' occurred at WCYN first; not a real event     
                continue
            
            time = event.get('time')
            if time[0].hour < 8 or time[0].hour > 18:
                # 'event' occurred before 8am or after 6pm; not a real event
                continue
            
            timeDiff = time[1] - time[0] # time difference in seconds          
            if timeDiff > (tpred+5) or timeDiff < (tpred-3):
                # 'event' did not reach the second station within predicted time
                #TODO: arbitrarially chose a window of 5/3s around predicted time
                continue
            
            # otherwise, this is a real event!
            realTrigDay.append(event)
            
            # make a 1 min window around event
            dt = datetime.timedelta(seconds=0,minutes=1)
            
            # save data in window around event
            event['data1'] = dataToSave(data1, time[0]-dt, time[0]+dt)
            event['data2'] = dataToSave(data2, time[1]-dt, time[1]+dt)
            
            # save event ID
            event['ID'] = ID
            ID += 1
        
        # append events from one day to dictionary
        trig[day] = realTrigDay
        print('Events today: ', len(realTrigDay))
        
    return trig


def main():
    # path to data (mseed)
    path1 = '/home/mad/Documents/Research2020/MSEEDdata/MMTN/'
    path2 = '/home/mad/Documents/Research2020/MSEEDdata/WCYN/'
    
    # path to save
    save = '/home/mad/Documents/Research2020/trig/'
    
    # [STA window length (sec), LTA window length (sec)]
    w_s = [2,30]
    
    # [lower, upper threshold values where trigger is deactivated, activated]
    thres = [1,5]
    
    dec = 5
    detrend = 'linear'
    filt = 'bandpass'
    
    # for filtering with bandpass
    freqmin = 1
    freqmax = 7
    
    # number of stations
    coSumThres = 2
    
    trig = pickEvents(path1,path2, w_s, thres, coSumThres, dec, detrend, filt, freqmin=freqmin,freqmax=freqmax)
    
    # save event triggers and data to load later
    np.save((save+'eventTrig_subset.npy'), trig) 
    
    return



main()



