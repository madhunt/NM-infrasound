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
#import time


# def func

path = '/home/mad/Documents/Research2020/MSEEDdata/MMTN/' # path to data (mseed)
w_s = [2,20] # [STA window length (sec), LTA window length (sec)]

thres = [0.5,5] # [low]

#thres_low = 0.5 # lower threshold value where trigger is deactivated below this point
#thres_high = 5 # upper threshold value where trigger is activated above this point

# for filtering with bandpass
freqmin = 1 # lower cutoff frequency
freqmax=10 # upper cutoff frequency

##TODO figure out how to make this automatic
npoints = 3600 * 250

##TODO: turn this whole file into a function
#def pickevents(path, w_sec, thres, freq):


# In[]: function to filter and detrend data, and compute STA/LTA ratios
def filt_stalta(data, w_s, detrend, filt, **filt_kwargs):
    '''
    Filters and detrends data, and computes STA/LTA ratios
    # INPUT:
        # data = trace
        # w_s = [STA,LTA] window sizes in sec
        # detrend = string; 'simple','linear','constant','polynomial','spline'
        # filt = type of filter; 'bandpass','lowpass','highpass', etc.
        # filt_kwargs = options depending on filter chosen (eg. freqmin,freqmax)
            # NOTE: df not needed as input
    # RETURNS: array of STALTA ratios from filtered datas
    '''
    
    # STA and LTA window lengths in terms of data points
    samprate = data.stats.sampling_rate
    STA_w = (samprate * w_s[0])
    LTA_w = (samprate * w_s[1])
    
    # center data at 0
    data_filt = data.copy()
    data_filt.detrend(detrend)
    # filter data
    data_filt.filter(filt,**filt_kwargs)
    
    #st=time.time()
    # calculate STA/LTA ratios
    cft = tg.classic_sta_lta(data_filt,STA_w,LTA_w)
    #print("classic_sta_lta_py call time: ",time.time()-st)
    return cft,samprate

# In[]: use filt_stalta on data

# make a list of data files in directory
files = [f for f in listdir(path) if isfile(join(path, f))]

# empty array to store events data
events = np.array([])

for filenum in range(0,np.size(files)+1):
    print('File Number '+str(filenum)+' Out of '+str(np.size(files)))
    
    # read in the data (in mseed format)
    mseed = obs.read(join(path, files[filenum]))
    data = mseed[0]
    
    cft,samprate = filt_stalta(data,w_s,'linear','bandpass',freqmin=freqmin,freqmax=freqmax)
    
    # if there is not a full hour of data, pad end with zeros
    if np.size(cft) < npoints:
        npad = npoints-np.size(cft)
        cft = np.pad(cft, npad, 'constant')
    
    ##TODO figure out why some files are too big
    if np.size(cft) > npoints:
        # there is some dumb poop in this file (1080250 points?????)
        print('uh oh... this file is too big???')
        continue
   
    assert np.size(cft)==samprate*3600, ('Sample rate may have changed...')
    
    # choose events
    new_events = (tg.trigger_onset(cft,thres[1],thres[0]), files[filenum])
    
    events = np.append(events, new_events)

# In[]: choose events and check

# choose events

t_st = int(files[0].split('.')[2]) # start time hour

# do a reality check
def getHour(event):
    t_hr = (t_st + event[0] / npoints) % 24
    return t_hr

for event in events:
    assert getHour(event) < 20, ('This event is after 8pm.. that is not right')







# In[]: plot

##TODO
# view filtered data with one day plot (ask oliver for plotting triggers code)
    # merge data





# In[]: TESTING ABOVE ON ONE FILE
# make a list of files
files = [f for f in listdir(path) if isfile(join(path, f))]

filenum=50

# empty array to store all cft data (each timeseries in a row)
all_cft = np.array([])

print('File Number '+str(filenum)+' Out of '+str(np.size(files)))

# read in the data (in mseed format)
mseed = obs.read(join(path, files[filenum]))
data = mseed[0]

cft = filt_stalta(data,w_s,'linear','bandpass',freqmin=freqmin,freqmax=freqmax)

# add data to all_cft matrix
all_cft = np.append(all_cft, cft)
    

# choose events
events = tg.trigger_onset(all_cft,thres[1],thres[0])



# In[]: check events (none at night, etc)








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
    
    ##TODO: ask about tapering
    data_filt.taper(max_percentage=0.05, type='hann')
    
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
    
    
    

# pull out ratio of amps -- window around each station
# no cammon at night or lunch -- every 15-30 min
# make a google earth map with the station locations on it
# compare with quarry blasts

# cluster analysis
