#!/usr/bin/env python3
# short term average long term average automatic event picker
# madhunt
# created 31 Jan 2020

import numpy as np
import matplotlib.pyplot as plt
import obspy as obs
import obspy.signal.trigger as tg
import obspy.signal.filter as fl
from os import listdir
from os.path import isfile, join
#import time # for debugging purposes

# In[]: define funtions
###############################################################################
def stalta(data, w_s, dec, detrend, filt, **filt_kwargs):
    '''
    Decimates, detrends, and filters data, and computes STA/LTA ratios.
    INPUTS:
        data: data in a trace
            type: trace
        w_s: [STA,LTA] window sizes in seconds
            type: 2 value array
        dec: factor to downsample the data by
            type: integer
        detrend: 'simple','linear','constant','polynomial','spline'
            type: string
        filt: type of filter; 'bandpass','lowpass','highpass'
            type: string
        filt_kwargs: options depending on filter chosen (eg. freqmin,freqmax)
            NOTE: df not needed as input
    RETURNS: array of STALTA ratios from processed data
    '''
    # STA and LTA window lengths in terms of data points
    samprate = data.stats.sampling_rate
    STA_w = (samprate * w_s[0])
    LTA_w = (samprate * w_s[1])
    
    # decimate data (decrease sampling rate)
    data_proc = data.copy()
    data_proc.decimate(dec)
    # center data at 0
    data_proc.detrend(detrend)
    # filter data
    data_proc.filter(filt,**filt_kwargs)
    
    # calculate STA/LTA ratios
    cft = tg.classic_sta_lta(data_proc,STA_w,LTA_w)
    
    return cft,samprate

###############################################################################
def pickEvents(path_data, path_save, w_s, thres, dec, detrend, filt, **filt_kwargs):
    '''
    Loads in data, runs stalta function, saves events in .npy file
    INPUTS:
        path_data: path to data (stored in mseed format)
            type: string
        path_save: path to save events (saved as .npy file)
            type: string
        w_s: [STA,LTA] window sizes in seconds
            type: 2 entry array
        thres: [lower, upper] threshold where trigger is deactivated or activated
            type: 2 entry array
        dec: factor to downsample the data by
            type: integer
        detrend: 'simple','linear','constant','polynomial','spline'
            type: string
        filt: type of filter; 'bandpass','lowpass','highpass', etc.
            type: string
        filt_kwargs: options depending on filter chosen
            if filt='bandpass': freqmin=lower cutoff frequency, freqmax=upper cutoff
            if filt='lowpass': freq=cutoff frequency
            if filt='highpass': freq=cutoff frequency
    RETURNS: 
    '''
    # make a list of data files in directory
    files = [f for f in listdir(path_data) if isfile(join(path_data, f))]
    
    events = np.array([]) # empty array to store events data
    
    # loop through each file
    for filenum in range(0,np.size(files)):
        print('File Number '+str(filenum)+' Out of '+str(np.size(files)))
        
        # read in the data
        mseed = obs.read(join(path_data, files[filenum]))
        data = mseed[0]
        
        cft,samprate = stalta(data, w_s, dec, detrend, filt, **filt_kwargs)
        
        npoints = 3600 * samprate / dec
        
        # if there is not a full hour of data, pad end with zeros
        if np.size(cft) < npoints:
            npad = npoints-np.size(cft)
            cft = np.pad(cft, npad, 'constant')
        
        ############################################
        ##TODO figure out why some files are too big
        if np.size(cft) > npoints:
            # there is something weird in this file (135 had 1080250 points?????)
            print('uh oh... file number ' + str(filenum) + 'is too big???')
            continue
        ############################################ 
        
        # make sure that the sample rate is the same and the array is the expected length
        assert np.size(cft)==samprate*3600, ('Sample rate may have changed...')
        
        # choose events
        new_events = (tg.trigger_onset(cft,thres[1],thres[0]))#, files[filenum])
        
        # events array 
        events = np.append(events, new_events)
    
    # save events
    np.save((path_save + 'events'), events)
    
    return events

###############################################################################
def checkEvents():
    
    ##TODO: finish checks
    
    t_st = int(files[0].split('.')[2]) # start time hour
    
    def getHour(event):
        '''Find start hour of event'''
        t_hr = (t_st + event / npoints)%24
        return t_hr
    
    # check that events are at reasonable times
    for event in events:
        print(event)
        assert getHour(event) < 20, ('This event is after 8pm.. that is not right')
        #assert getHour(event) != 12, ('This event occured during lunch.. that is strange')
    
    # compare events between stations
    
    
    
    return

###############################################################################
def plotEvents():
    
    ##TODO
    
    # save plots
    
    return



# In[]: call functions
path_data = '/home/mad/Documents/Research2020/MSEEDdata/MMTN/' # path to data (mseed)
path_save = '/home/mad/Documents/Research2020/picks/MMTN' # no final / so MMTN is in filename

w_s = [2,20] # [STA window length (sec), LTA window length (sec)]

thres = [0.5,5] # [lower, upper threshold values where trigger is deactivated, activated]

dec = 5
detrend = 'linear'
filt = 'bandpass'

# for filtering with bandpass
freqmin = 1
freqmax = 7

# pick events
events = pickEvents(path_data, path_save, w_s, thres, dec, detrend, filt, freqmin,freqmax)



# In[]: compare events with amrit's picks




# In[]: one day plot of events

#files = [f for f in listdir(path_data) if isfile(join(path_data, f))]
#mseed = obs.read(join(path_data, files[0]))
#
#for filenum in range(1,24): # read only files for the first day
#    
#    mseed += obs.read(join(path_data, files[filenum]))
#
#mseed.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
#mseed.plot(type='dayplot')



# In[]: example for one hour
def ex_1hr():
    STA_ws = 5 # short term window length (s)
    LTA_ws = 50 # long term window length (s)
        
    thres_high = 5 # upper threshold value where trigger is activated above this point
    thres_low = 0.5 # lower threshold value where trigger is deactivated below this point
    
    # make a list of files
    files = [f for f in listdir(path) if isfile(join(path, f))]
    
    filenum=50
    mseed = obs.read(join(path, files[filenum]))
    data = mseed[0]
    amp = data.data
    samprate = data.stats.sampling_rate
    
    # STA and LTA window lengths in terms of data points
    STA_w = samprate * STA_ws
    LTA_w = samprate * LTA_ws
    
    # center data at 0
    data_filt = data.copy()
    data_filt.detrend('linear')
    
    # filter data
    data_filt = fl.bandpass(data_filt,freqmin=1,freqmax=10,df=250)
    
    # take sta/lta ratio
    cft = tg.classic_sta_lta(data_filt,STA_w,LTA_w)
    
    # find trigger on/off times given the higher/lower thresholds
    on_off = tg.trigger_onset(cft, thres_high, thres_low)
    
    #plt.figure(1)
    ## plotting code modified from https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html
    #ax = plt.subplot(211)
    #plt.title('Filtered Signal')
    #plt.plot(data_filt, 'k')
    #ymin, ymax = ax.get_ylim()
    #plt.vlines(on_off[:, 0], ymin, ymax, color='r', linewidth=2) # plot trigger on
    #plt.vlines(on_off[:, 1], ymin, ymax, color='b', linewidth=2) # plot trigger off
    #plt.subplot(212, sharex=ax)
    #plt.title('STA/LTA Ratio')
    #plt.plot(cft, 'k')
    #plt.hlines([thres_high, thres_low], 0, len(cft), color=['r', 'b'], linestyle='--')
    #plt.axis('tight')
    #plt.suptitle(files[filenum]) # make a title with the filename
    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.show()
    
    
    plt.figure(1)
    # plotting code modified from https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html
    ax = plt.subplot(311)
    plt.title('Filtered Signal')
    plt.plot(data_filt, 'k')
    ymin, ymax = ax.get_ylim()
    plt.vlines(on_off[:, 0], ymin, ymax, color='r', linewidth=2) # plot trigger on
    plt.vlines(on_off[:, 1], ymin, ymax, color='b', linewidth=2) # plot trigger off
    plt.subplot(312, sharex=ax)
    plt.title('STA/LTA Ratio')
    plt.plot(cft, 'k')
    plt.hlines([thres_high, thres_low], 0, len(cft), color=['r', 'b'], linestyle='--')
    
    ax = plt.subplot(313, sharex=ax)
    plt.title('Unfiltered Signal')
    plt.plot(data, 'k')
    
    plt.axis('tight')
    plt.suptitle(files[filenum]) # make a title with the filename
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    