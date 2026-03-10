import numpy as np

def square_wave(tau_mag, t, beg, end):
    tau = np.zeros(len(t))
    tau = np.array(tau)
    start = np.argmin(np.fabs(t - beg))
    end = np.argmin(np.fabs(t - end))
    for i in range(start,end):
        tau[i] = tau_mag

    return tau

def error(fit, solution, flag):
    if(flag == 0):
        err = np.abs(solution - fit)/np.abs(solution)
    else:
        err = np.linalg.norm(fit - solution, 2) / np.max(solution)

    return err

    
