import numpy as np

def combine_data(full_t, full_y, index, phi_an, proc_data_switch, N_f, N_f2):
    half = N_f2
    
    if(proc_data_switch == 1):
        full_y = full_y - full_y[0]  # clean out some data

        spot = full_t[index]  # t value at peak
    

    dt = np.mean(np.diff(full_t))

    if(proc_data_switch == 1) or (len(full_t) < half):
        # define end of positive t
        n_end = half - len(full_t[:index])
        t_end = np.linspace(spot + dt, spot + dt * n_end, n_end-1)
    else:
        t_end = None
    
    if(proc_data_switch == 1) or (len(full_t) < half):
        # define negative t
        t_neg = np.linspace(-half * dt, -dt, half)
        y_neg = np.zeros(len(t_neg))
    else:
        t_neg = np.linspace(-half * dt, -dt, half)
        y_neg = np.full(half, full_y[0])


    #print(len(t_neg), len(full_t[:index]), len(t_end))
    #print(len(y_neg), len(full_y[:index]), len(phi_an(t_end)))

    if(proc_data_switch == 1):
        # combine
        franken_t = np.concatenate([t_neg, full_t[:index+1], t_end])
        franken_y = np.concatenate([y_neg, full_y[:index+1], phi_an(t_end)])
    elif(len(full_t) < half):
        tail = np.zeros(len(t_end))
        franken_t = np.concatenate([t_neg, full_t[:index+1], t_end])
        franken_y = np.concatenate([y_neg, full_y[:index+1], tail])
    else:
        franken_t = np.concatenate([t_neg, full_t[:half]])
        franken_y = np.concatenate([y_neg, full_y[:half]])

    return franken_t, franken_y, t_end


def remove_noise(full_t, full_y, threshold):
    
    dy = np.abs(np.diff(full_y) / np.diff(full_t))  # derivative
    start_index = np.argmax(dy > threshold)

    start_mean = np.mean(full_y[:start_index])
    full_y -= start_mean
    
    for i in range(start_index):
        full_y[i] = 0

    return full_y


import numpy as np
from scipy.signal import butter, filtfilt, welch

def denoise_from_quiet_region(data, sample_rate, quiet_end_idx, pad_zeros=True):
    """
    Characterizes noise from a quiet (should-be-zero) region at the start,
    then filters the entire signal to remove that noise.
    
    Parameters
    ----------
    data         : 1D np.ndarray
    sample_rate  : float, samples per second (or 1/dt)
    quiet_end_idx: int, last index of the quiet region
    pad_zeros    : bool, if True pads zeros at both ends before filtering
    
    Returns
    -------
    cleaned : 1D np.ndarray, same length as data
    cutoff  : float, the noise cutoff frequency that was found
    """
    quiet = data[:quiet_end_idx]

    # Estimate noise power spectrum in the quiet region
    freqs, power = welch(quiet, fs=sample_rate, nperseg=min(len(quiet), 64))

    # Find the dominant noise frequency (peak power in quiet region)
    dominant_noise_freq = freqs[np.argmax(power)]

    # Set cutoff just below the noise — use 80% of dominant freq as cutoff
    # (keeps signal content below it, kills noise at and above it)
    cutoff = dominant_noise_freq * 0.3
    cutoff = max(cutoff, freqs[1])  # guard against cutoff=0

    print(f"Dominant noise frequency: {dominant_noise_freq:.4f}")
    print(f"Cutoff set to:            {cutoff:.4f}")

    # Design zero-phase Butterworth low-pass filter
    nyq = 0.5 * sample_rate
    order = 4  # good default: steep enough without ringing
    b, a = butter(order, cutoff / nyq, btype='low')

    # Pad with zeros at both ends to avoid edge artifacts
    if pad_zeros:
        pad = len(data) // 4  # pad with 25% of signal length
        padded = np.concatenate([np.zeros(pad), data, np.zeros(pad)])
        filtered = filtfilt(b, a, padded)
        cleaned = filtered[pad:-pad]
    else:
        cleaned = filtfilt(b, a, data)

    return cleaned, cutoff
