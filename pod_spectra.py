import numpy as np
import time
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats

import sys
sys.path.insert(0, "/project/def-alamj/shared/include/pyutils")

import awcmviewerutils as av

"""
data collection tool
/project/def-alamj/shared/bin/v2306/xpostabl -INP post_probes.inp -analysis probes -fields 'U' -time 500 -out fowf15mwR0 -pwd -pod 1 # raw data
/project/def-alamj/shared/bin/v2306/xpostabl -INP post_probes.inp -analysis probes -fields 'U' -time 500 -out fowf15mwR0 -pwd -pod 0 # rank 1 pod
/project/def-alamj/shared/bin/v2306/xpostabl -INP post_probes.inp -analysis probes -fields 'U' -time 500 -out fowf15mwR0 -pwd -pod -1 # rank 1 pod residual

"""


"""
adjust the user dependent path.

fowf15mwR2prb.h5 - original time series, u
fowf15mwR2p1b.h5 - POD filtered time series, u - u(pod)
fowf15mwR2p0b.h5 - POD filtered time series, u(pod)

This example code illustrates effects of POD on the spectra and autocorrelation
"""

ftag = "fowf15mwR2"
base = "/scratch/alamj/WindFarms/"
out1 = base + "fowf15mwR2/statistics/xdmf/fowf15mwR2prb.h5" # raw data
out2 = base + "fowf15mwR2/statistics/xdmf/fowf15mwR2p1b.h5" # rank1 residual
out0 = base + "fowf15mwR2/statistics/xdmf/fowf15mwR2p0b.h5" # rank1 reconstruction
print("POD rank 1 file: ", out0)
print("Original signals file: ", out1)
print("POD rank 1 residual file: ", out2)

xdmf = base + "fowf15mwR2/statistics/xdmf/"
jpg = xdmf + "fowf15mwR2Ek.jpg"

# target probe location according your configuration

#loc = [10, 3720, 10]
loc = [3360, 1320, 150]
#loc = [3360, 2520, 150]
#loc = [3360, 3720, 150]
#loc = [3360, 4920, 150]

# original time series
# read time
t1 = av.fetch_time_series(out1, loc, 0, "time")
# read velocity
u1 = av.fetch_time_series(out1, loc, 0, "U")
# read spectrum
Ek1= av.fetch_time_series(out1, loc, 0, "psd")
# read frequency
fq1 = av.fetch_time_series(out1, loc, 0, "freq")

# read u-u(pod), rank-1 pod residual
# read time
t2 = av.fetch_time_series(out2, loc, 0, "time")
# read velocity
u2 = av.fetch_time_series(out2, loc, 0, "U")
# read spectrum
Ek2= av.fetch_time_series(out2, loc, 0, "psd")
# read frequency
fq2 = av.fetch_time_series(out2, loc, 0, "freq")

# read rank1 POD reconstruction
# read spectrum
Ek0 = av.fetch_time_series(out0, loc, 0, "psd")

speed = 10
# time step, assume the same time step for both data, if not, adjust accordingly
dt = t1[1]-t1[0]
n = len(t1) # number of time steps

df = 1/(len(Ek1)*dt) # frequency step
fq = fq1*df          # convert sampling frequency
fq = 2.0*np.pi*fq/speed # normalize
k = fq1


plt.figure()
# plot the orginal spectra
plt.loglog(fq, Ek1/(n*dt), 'b-', label='orginal')
# plot the rank-1 POD residual spectra
plt.loglog(fq, Ek2/(n*dt), 'k--', label='POD-1')

plt.loglog(fq[10:1000], 1e5*k[10:1000]**(-5/3),'ro',label='-5/3')
plt.loglog(fq[10:], 5e1*k[10:]**(-5/3),'r-',label='-5/3')

"""
In the simulation, dx=dy=dz = 8 m.
So, theoretical cutoff is pi/8.
Due to test filter at 2*dx, expected cutoff is pi/(2*8)
If the diamond indicates the drop, it suggest expected filtering accuracy
"""
plt.loglog(0.2,10,'kd',label='cutoff k=0.2')
plt.gca().set_box_aspect(1)
plt.ylim([1e-10,1e10])
plt.xlabel('k')
plt.ylabel('Ek')
plt.legend()

print(jpg)
plt.savefig(jpg, bbox_inches='tight')
plt.close()


"""
In the following code, we compare the time series of
original u
residual v = u - u_reconstructed
"""


# detrend (find linear trend)
# find a linear curve fit of u as m*t + c
A = np.vstack([t1 - t1.mean(), np.ones_like(t1)]).T # n x 2 matrix

# solve a least square problem to evaluate slope m, y-intercept c
m, c = np.linalg.lstsq(A, u1, rcond=None)[0]
#u_d = u - (a*(t - t.mean()) + b)  # detrended
# best linear fit
u_l = (m*(t1 - t1.mean()) + c)


plt.figure()
# for simplicity, use t = t1
t = t1
# original time series
plt.plot(t-t[0], u1, 'b--', label='original')
# residual of rank-1 POD
plt.plot(t-t[0], u2, 'r--', label='POD-1')
# rank-1 POD
plt.plot(t-t[0], u1-u2, 'k-', label='POD trend')
plt.plot(t-t[0], u_l, 'g--', label='best linear fit')

plt.legend()
plt.savefig(xdmf + "fowf15mwR2ts.jpg", bbox_inches='tight')
plt.close()

"""
Following code computes autocorrelation

"""
# tau1, R1 = av.get_acf(u1, dt) # test using time series
tau2, R2 = av.psd_to_acf(fq, Ek2, dt)
tau0, R0 = av.psd_to_acf(fq, Ek0, dt)

# auto-correlation is done

Uc = np.mean(u1-u2)
print("Mean = ", np.mean(u1), "POD mean = ", Uc)

print("\nStrongest POD mode stat (large scale):")
res = av.get_coherent_scales(tau0, R0, Uc)
print("Dominant period = ", round(res['T_p']+0.5), "[s], Coherent time = ",  round(res['T_c']+0.5), "[s], Dominant length = ", round(res['Lambda_p']+0.5), "[m], Coherent length = ", round(res['Lambda_c']+0.5), "[m]")

print("\nStrongest POD residual stat (small scale):")
res = av.get_coherent_scales(tau2, R2, Uc)
print("Dominant period = ", round(res['T_p']+0.5), "[s], Coherent time = ",  round(res['T_c']+0.5), "[s], Dominant length = ", round(res['Lambda_p']+0.5), "[m], Coherent length = ", round(res['Lambda_c']+0.5), "[m]")

plt.figure()
plt.plot(tau0, R0, 'k--', label='Ruu, rank 1 POD trend')
plt.plot(tau2, R2, 'b-', label='Ruu, rank 1 POD residual')

plt.legend()

plt.savefig(xdmf + "fowf15mwR2acf.jpg", bbox_inches='tight')
plt.close()

"""
Axial induction and thrust and flap-wise bending moment

Here, we estimate an equivalent per-blade flapwise root moment
from actuator-disk velocities via momentum theory.

Using a(t) 1 ud(t)/speed,
thrust, T(t) = 2\rho Aud^2(t) a(t)/(1-a(t)).

Then, M_flap = 3/(3B)Rd/2

"""

Rd = 240
CT = 0.560001
CP = 0.465732;
a_ind = (1-CP/CT)
ud = speed*(1-a_ind)
# rated power
w_pwr = av.compute_power(ud, 240, CP, CT, 1.23)/1e6
# time series power by LES
pwr = av.compute_power(u1-u2, 240, CP, CT, 1.23)/1e6

#axial induction time series
a = 1-(u1-u2)/speed
print("induction ", a_ind, "mean induction ", np.mean(a), "\nrated power = ", w_pwr, np.mean(pwr))

T1 = 2.0*1.23*np.pi*((Rd/2)**2)*(u1)**2*a/(1-a)
T = 2.0*1.23*np.pi*((Rd/2)**2)*(u1-u2)**2*a/(1-a)
Mflap1 = 2.0/(3*3) * T1 * Rd/2
Mflap = 2.0/(3*3) * T * Rd/2
plt.figure()
plt.plot(t-t[0], Mflap, 'b--')
plt.plot(t-t[0], Mflap1, 'k--')

plt.savefig(xdmf + ftag + "axl.jpg", bbox_inches='tight')
plt.close()



"""
The following code computes the rank-1 POD using python code.

"""
exit()

print('Compute POD: ', out1)
Vr = av.get_by_pod_rank(out1, "U", 1)

# select the index of desired data at probe loc
id = av.get_row_id(out1, loc, 3, 0) # 3 component field, look for 0-th component

plt.figure()
# plot POD reconstruction by python
plt.plot(t1-t1[0], Vr[id,:], 'k--', label='python POD')

# plot POD reconstruction by C++
plt.plot(t1-t1[0], u1-u2, 'r--', label='c++ POD')

plt.legend()
plt.savefig(xdmf + "fowf15mwR2tsRec.jpg", bbox_inches='tight')
plt.close()
