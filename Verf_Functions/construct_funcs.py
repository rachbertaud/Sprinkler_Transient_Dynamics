import numpy as np

def square_wave(tau_mag, t, beg, end):
    tau = np.zeros((len(t), 1))
    start = np.argmin(np.fabs(t - beg))
    end = np.argmin(np.fabs(t - end))
    for i in range(start,end):
        tau[i] = tau_mag

    return tau
