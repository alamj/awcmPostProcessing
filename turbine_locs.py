import numpy as np
import h5py
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "/project/def-alamj/shared/include/pyutils")

import awcmviewerutils as av


loc = av.get_turbine_location(3360, 1320)

val = np.random.randint(5, 15, size=loc.shape[0])

for i in range(loc.shape[0]):
    p = loc[i]
    if p[2] == 90:
        val[i] = 5
    else:
        val[i] = 15
    
av.plot_turbine(loc,val)
