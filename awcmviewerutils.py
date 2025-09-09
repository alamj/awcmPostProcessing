import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import correlate


# read an attribute from .h5
def get_data(hdf, attr):

    f = h5py.File(hdf, 'r')
    # Check if it's a dataset
    if attr in f:
        return f[attr][:]
    # Check if it's an attribute in root group
    elif attr in f.attrs:
        return f.attrs[attr]
    else:
        return None


def fetch_time_series(probe, probe_point, comp, data="time"):
    '''
    probe = /path/to/simulation/statistics/xdmf/filename.h5
    probe_point = [x,y,z], which maybe very close to the probe point, 
    ideally, (x,y,z) should be exactly the desired probe location
    n_comp = 1 for scalar 3 for vector 6 for symmetric tensor
    comp is component of the data, i.e. 0 for Ux and so on
    '''

    n_comp = 3
    if data == "UPrime2Mean":
        n_comp = 6
        
    if data == "time":
        t = get_data(probe, data).reshape(-1)
        return t
    
    points = get_data(probe, "points")
    ptr = np.array(probe_point)

    u = get_data(probe, data)

    if u is None:
        print("Incorrect query: ", data)
        return None
    
    distances = np.linalg.norm(points - ptr, axis=1)
    idx = np.argmin(distances)
    print("probe id: ", idx, " probe location: ", points[idx])

    id = (idx)*n_comp + comp
    print("row index ", id)
    u = u[id,:]

    return u


def get_time_series(probe, probe_point, n_comp, comp, data="data"):
    '''
    probe = /path/to/simulation/statistics/xdmf/filename.h5
    probe_point = [x,y,z], which maybe very close to the probe point, 
    ideally, (x,y,z) should be exactly the desired probe location
    n_comp = 1 for scalar 3 for vector 6 for symmetric tensor
    comp is component of the data, i.e. 0 for Ux and so on
    '''

    if comp > n_comp -1 :
        print ("Error: comp ", comp, "should be reduced")
        return 0
    
    points = get_data(probe, "points")
    ptr = np.array(probe_point)

    u = get_data(probe, data)

    if u is None:
        print("Incorrect query: ", data)
        return None, None
    
    t = get_data(probe, "time")
    if t is None:
        t = u[0,:]
        u = u[1:,:]
    else:
        t = t.reshape(-1)

    distances = np.linalg.norm(points - ptr, axis=1)
    idx = np.argmin(distances)
    print("probe id: ", idx, " probe location: ", points[idx])

    id = (idx)*n_comp + comp
    print("row index ", id)
    u = u[id,:]

    return t, u 

def get_row_id(probe, probe_point, n_comp, comp):
    points = get_data(probe, "points")
    ptr = np.array(probe_point)
    
    distances = np.linalg.norm(points - ptr, axis=1)
    idx = np.argmin(distances)
    

    id = (idx)*n_comp + comp

    return id


# compute energy spectrum
def get_spectra(u, dt=1.0):
    """
    Compute PSD of each row in u using numpy.fft.fft
    Parameters:
        u  : 2D numpy array (m x n), where each row is a time series of length n
        dt : time step (default = 1.0)
    Returns:
        out : 2D numpy array of shape (m+1, n//2 + 1)
              out[0] is frequency array
              out[1:] are PSDs of each row in u
    """
    if u.ndim == 1:
        u = u.reshape((1,u.shape[0]))
        
    m, n = u.shape
    freqs = np.fft.rfftfreq(n, d=dt)
    out = np.zeros((m + 1, len(freqs)))

    out[0] = freqs

    for i in range(m):
        fft_vals = np.fft.rfft(u[i])
        psd = (np.abs(fft_vals) ** 2) / (n * dt)
        out[i + 1] = psd

    return out


# two sided ccf
def awcm_ccf_two_sided(Ulow, Uhgh, fs=1):
    ccf_two_sided = correlate(Ulow - Ulow.mean(), Uhgh - Uhgh.mean(), mode='full', method='auto')
    lags = np.arange(-len(Ulow) + 1, len(Ulow)) / fs  # fs = sampling frequency
    ccf_two_sided /= (np.std(Ulow) * np.std(Uhgh) * len(Ulow))  # normalize
    r_peak = np.argmax(ccf_two_sided)
    l_peak = np.argmin(ccf_two_sided)
    print('left peak ', l_peak, ', right peak', r_peak)

    return lags, ccf_two_sided, l_peak, r_peak


def sure_threshold(x, sigma):
    """
    Compute SURE-based soft threshold for a given detail coefficient vector x and noise sigma.
    """
    x = np.asarray(x)
    n = x.size
    x2 = np.sort(np.abs(x))**2
    risks = (n - 2*np.arange(n) + np.cumsum(x2)) / n
    idx_min = np.argmin(risks)
    threshold = np.sqrt(x2[idx_min])
    return threshold

def wavelet_sure_filter(signal, wavelet='sym8', level=3):
    """
    Apply SURE-based soft thresholding wavelet denoising to a 1D signal.

    Parameters:
    -----------
    signal : ndarray
        1D input signal (e.g., wind velocity).
    wavelet : str
        Wavelet type (default: 'sym8').
    level : int
        Decomposition level (default: 3).

    Returns:
    --------
    signal_denoised : ndarray
        Wavelet-filtered (denoised) signal.
    """
    signal = np.asarray(signal).flatten()
    sigma = estimate_sigma(signal, average_sigmas=True)
    coeffs = pywt.wavedec(signal, wavelet=wavelet, level=level)

    coeffs_filtered = [coeffs[0]]  # Keep approximation (low-freq) unchanged
    #for i, cD in enumerate(coeffs[1:], start=1):
        #threshold = sure_threshold(cD, sigma)
    
        # Apply extra attenuation to low-frequency levels
        # Example: scale threshold for coarser levels
        #scaling_factor = 1.0 + 0.5 * (level - i)  # strongest for low levels
        #threshold *= scaling_factor

    for cD in coeffs[1:]:
        threshold = sure_threshold(cD, sigma)
        cD_thresh = pywt.threshold(cD, threshold, mode='soft')
        coeffs_filtered.append(cD_thresh)

    signal_denoised = pywt.waverec(coeffs_filtered, wavelet=wavelet)
    return signal_denoised[:len(signal)]




def bandpass(u, fs, lowcut=0.05, highcut=0.5, order=4):
    b, a = butter(order, [lowcut / (0.5 * fs), highcut / (0.5 * fs)], btype='band')
    return filtfilt(b, a, u)



def get_turbine_location(startX, startY, hub1=150, hub2=90, D=240, n_pairs=5, spanwise_spacing=5, streamwise_spacing=4):
    """
    Returns turbine locations in a staggered layout with 5 extra turbines at the end (hub1 height).
    Even x-columns: 5 turbines at hub1 (150m)
    Odd x-columns:  4 staggered turbines at hub2 (90m)
    """
    dx = streamwise_spacing * D
    dy = spanwise_spacing * D

    locs = []
    for j in range(n_pairs * 2):  # 10 columns: 5 even, 5 odd
        x = startX + j * dx
        if j % 2 == 0:
            # 5 turbines at hub1
            for i in range(5):
                y = startY + i * dy
                locs.append([x, y, hub1])
        else:
            # 4 staggered turbines at hub2
            for i in range(4):
                y = startY + dy / 2 + i * dy
                locs.append([x, y, hub2])

    # Add 5 extra turbines at the end, fixed x, increasing y
    x_extra = startX + (2 * n_pairs) * dx  # Continue x progression
    for i in range(5):
        y = startY + i * dy
        locs.append([x_extra, y, hub1])

    return np.array(locs)


def plot_turbine(locs, vals):
    """
    Plots turbine layout based on 3-column locs: [x, y, hub_height]
    Red circles for hub = 150 m, blue stars for hub = 90 m
    Labels each turbine with the corresponding entry in `vals`.
    """
    locs = np.array(locs)

    # Calculate plot limits with 10% margin
    x_vals = locs[:, 0]
    y_vals = locs[:, 1]

    x_margin = 0.2 * (x_vals.max() - x_vals.min())
    y_margin = 0.4 * (y_vals.max() - y_vals.min())

    xmin = x_vals.min() - x_margin
    xmax = x_vals.max() + x_margin
    ymin = y_vals.min() - y_margin
    ymax = y_vals.max() + y_margin

    # Separate by hub height
    locs_150 = locs[locs[:, 2] == 150]
    locs_90 = locs[locs[:, 2] == 90]

    plt.figure(figsize=(10, 6))
    plt.plot(locs_150[:, 0], locs_150[:, 1], 'ro', label='Hub 150 m')
    plt.plot(locs_90[:, 0], locs_90[:, 1], 'b*', label='Hub 90 m')

    for i, (x, y, _) in enumerate(locs):
        plt.text(x + 20, y, str(vals[i])+ ' MW', fontsize=10)

    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Turbine Layout (Red: 150 m, Blue: 90 m)")
    plt.grid(False)
    plt.legend()
    plt.axis('equal')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.tight_layout()
    plt.show()



def compute_power(u_mag, D=240, Cp=0.486, Ct=0.86, rho=1.225):
    """
    Compute instantaneous power P(t) from wind velocity time series u_t.

    Parameters:
    - u_mag(t): (N x 1) array of velocity components u at the rotor center
    - D: rotor diameter in meters
    - Cp: power coefficient
    - rho: air density (kg/m^3)

    Returns:
    - P(t): instantaneous power in Watts
    """

    a = 1 - Cp/Ct
    
    A = (np.pi * D**2) / 4  # Rotor swept area
    #u_mag[:] = np.mean(u_mag) # debug only remove later
    P = 2.0 * rho * A * u_mag**3 * a/(1-a)  # Instantaneous power

    return P



def epfl_data(exp, t, j):
    print('EPFL_EXP_location_'+str((t-1) * 5 + 1 + j)+'.txt')
    xEXP, yEXP = np.loadtxt(exp + 'EPFL_EXP_location_'+str((t-1) * 5 + 1 + j)+'.txt',skiprows=1, usecols=(0, 1), unpack=True)
    
    return xEXP, yEXP

