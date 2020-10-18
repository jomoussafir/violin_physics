import numpy as np

def smooth(x,h):
    """
    Basic symetric smoothing using uniform weights.
    
    w: numpy 1d array containing signal to be smoothed
    h: window length 
    """
    w = np.ones(h)
    w = w/sum(w)
    return np.convolve(x, w, mode="same")


def exp_smooth(x, h):
    phi = np.exp(np.log(0.5)/h)
    xs = np.zeros(x.size)
    
    xs[0]=x[0]
    for i in np.arange(1, x.size):
        xs[i] = phi*xs[i-1] + (1-phi)*x[i]

    xs_init = (xs - x[0]*phi**(x.size+1))/(1 - phi**(x.size+1))    
    return xs_init


def bnd(idx):
    d={}
    d["start"] = []
    d["end"] = []

    if idx[0]==1:
        d["start"].append(0)
    
    for i in range(1,len(idx)):
        if idx[i-1]==0 and idx[i]==1:
            d["start"].append(i)
        if idx[i-1]==1 and idx[i]==0:
            d["end"].append(i-1)

    if idx[-1]==1:
        d["end"].append(i)

    d["nb_interval"] = len(d["start"])
    return d


def wave_cut(wave, h, thr_ratio):
    """
    wave_cut makes a dictionnary of start and end times in seconds 
    from wave for removing silences
    
    h: smoothing parameter in seconds
    thr_ratio: if smoothed sound amplitude in dB goes below 
                thr_ratio*max amplitude, wave is considered as silence
    """
    
    ws = smooth(np.abs(wave.ys), int(h*wave.framerate))
    emin = np.quantile(np.log(ws), 0.05)
    emax = np.quantile(np.log(ws), 0.95)

    thr=emin + thr_ratio*(emax - emin)
    idx = np.array(np.log(ws) >= thr, dtype=int)
    t_dict=bnd(idx)
    
    for i in range(len(t_dict["start"])):
        t_dict["start"][i] = t_dict["start"][i]/wave.framerate

    for i in range(len(t_dict["end"])):
        t_dict["end"][i] = t_dict["end"][i]/wave.framerate

    return t_dict


def violin_spec(wave, fmin=20, fmax=5000):
        
    ## cut wave into small samples
    t_dict = wave_cut(wave, 0.5, 0.75)
    nb_sample = t_dict["nb_interval"]
    
    ## make segments
    segment = {}
    for i in range(t_dict["nb_interval"]):
        segment[i] = wave.segment(t_dict["start"][i], t_dict["end"][i] - t_dict["start"][i])
    
    ## compute spectrum for each segment
    spectrum = {}
    fs_sample = np.linspace(fmin, fmax, num=100000, endpoint=True)
    nb_fs_sample = len(fs_sample)

    ## make linear interpolation of spectrum for each spectrum
    ## with given frequencies
    spectrum_approx = np.zeros([nb_fs_sample,nb_sample])
    for i in range(nb_sample):
        spectrum[i] = segment[i].make_spectrum()
        x = spectrum[i].fs
        y = spectrum[i].amps
        f = interp1d(x,y)
        spectrum_approx[:,i] = f(fs_sample)
    
    ## compute average
    spectrum_approx_avg = np.mean(spectrum_approx, axis=1)
    
    ## make dict with average spectrum
    spec = {}
    spec["fs"] = fs_sample
    spec["amps"] = spectrum_approx_avg
    
    return spec



def peak(spec, thr=-0.1):
    """
    Compute location and amplitude of a spectrum peaks
    
    spec: a dictionnary with keys "fs" and "amps"
    
    Returns numpy array with 
    """
    x = spec["fs"]
    y = np.log(spec["amps"])
    Y = spec["amps"]
    
    ## zscore
    z = 2*(y-np.min(y))/(np.max(y) - np.min(y)) - 1

    zs = smooth(np.array(z>=thr, dtype=int),1000)
    d = bnd(zs>=0.05)
    
    nb_peak = d["nb_interval"]

    amax = np.zeros(nb_peak)
    fmax = np.zeros(nb_peak)
        
    for i in np.arange(nb_peak):
        amax[i] = np.max(Y[d["start"][i]:d["end"][i]])
        idx = np.where(Y[d["start"][i]:d["end"][i]] == amax[i])[0][0]
        fmax[i] = np.median(x[d["start"][i]:d["end"][i]][idx])

    return fmax, amax