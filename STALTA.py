# short term average through long term average automatic picker
# madhunt
# created 31 Jan 2020

# In[]: set up + define variables
import numpy as np
import obspy as obs
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

##TODO make this better (maybe read in all or one at a time but u shouldnt have to change it)
filenum = 80 # choose which of the files to read with this number
path = '/home/mad/Documents/Research2020/MSEEDdata/MMTN/' # path to data

threshold = 4.0 # threshold value of STA/LTA ratio
STA_del = 720 # short term window length
LTA_del = 72000 # long term window length

# In[]: load in the data (in mseed format)
files = [f for f in listdir(path) if isfile(join(path, f))]
mseed = obs.read(join(path, files[filenum]))
data = mseed[0].data

##TODO import the correct time sequence
l = 250 # currently using this to find time

# In[]: take moving averages
# takes a moving average given a 1D array of data and a window length
# returns an array of the calculated averages
def mov_avg(A, n):
    # A = array
    # n = window length
    C = np.cumsum(A)
    C[n:] = C[n:] - C[:-n]
    avgs = C[n-1:]/n
    return avgs

# take short term average
STA = mov_avg(data, STA_del)
t_STA = np.arange(0,l,l/np.size(STA))

# take long term average
LTA = mov_avg(data, LTA_del)
t_LTA = np.arange(0,l,l/np.size(LTA))


# In[]: take ratio
l = np.size(data)

l_STA = np.size(STA) #15
l_LTA = np.size(LTA) #4
l_rat = l_STA/l_LTA #15/4=3.75

v = np.arange(0,l_STA,l_rat)
print(v)



ratio = STA[1]/LTA[1]





# In[]: plotting yay
# now plot to see what it looks like
plt.plot(t_STA, STA, 'b')
plt.plot(t_LTA, LTA, 'r')
plt.show()




