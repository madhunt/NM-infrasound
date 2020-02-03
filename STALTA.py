# short term average through long term average automatic picker
# madhunt
# created 31 Jan 2020

# In[]: set up + define variables
import numpy as np
import obspy as obs
import obspy.signal.trigger as tg
#import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join


path = '/home/mad/Documents/Research2020/MSEEDdata/MMTN/' # path to data
STA_ws = 2 # short term window length (s)
LTA_ws = 10 # long term window length (s)


# In[]:

# make a list of files
files = [f for f in listdir(path) if isfile(join(path, f))]

# empty array to store all cft data (in each row)
all_cft = np.zeros((np.size(files),900000))


for filenum in range(0,np.size(files)):
    # read in the data (in mseed format)
    mseed = obs.read(join(path, files[filenum]))
    data = mseed[0]
    amp = data.data
    samprate = int(data.stats.sampling_rate)
    
    # get start and end time, convert from UTC to Unix time
    #t0 = int((data.stats.starttime).strftime('%s'))
    #tf = int((data.stats.endtime).strftime('%s'))
    

    # STA and LTA window lengths in terms of data points
    STA_w = samprate * STA_ws
    LTA_w = samprate * LTA_ws
    
    cft = tg.classic_sta_lta(amp,STA_w,LTA_w)
    tg.plot_trigger(data, cft, 1.5, 0.5)
    
    all_cft[filenum,:] = cft
    


# In[]: take moving averages
## takes a moving average given a 1D array of data and a window length
## returns an array of the calculated averages
#def mov_avg(A, n):
#    # A = array
#    # n = window length
#    C = np.cumsum(A)
#    C[n:] = C[n:] - C[:-n]
#    avgs = C[n-1:]/n
#    return avgs
#
## take short term average
#STA = mov_avg(amp, STA_w)
#t_STA = np.linspace(t0,tf,num=np.size(STA))
#
## take long term average
#LTA = mov_avg(amp, LTA_w)
#t_LTA = np.linspace(t0,tf,num=np.size(LTA))
#





