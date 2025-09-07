import numpy as np 
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats

import sys
sys.path.insert(0, "/project/def-alamj/shared/include/pyutils")

import awcmviewerutils as av


out = "/scratch/alamj/WindFarms/fowf15mwR0/statistics/xdmf/fowf15mwR0.h5"


t, u = av.get_time_series(out, [1680, 1320, 150], 3, 0, "U")

print(t.shape, u.shape)
plt.plot(t,u)
plt.show()}