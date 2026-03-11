import numpy as np

seed_value = 1010
rng = np.random.default_rng(seed_value)

def square_wave(tau_mag, t, beg, end, flag):
    tau = np.zeros(len(t))
    tau = np.array(tau)
    start = np.argmin(np.fabs(t - beg))
    end = np.argmin(np.fabs(t - end))
    len_t = end - start
    for i in range(start,end):
        tau[i] = tau_mag
    if(flag == 1):
        tau[start:end] += rng.normal(0,0.1*tau_mag,len_t)
        starter = int(np.floor(0.1*len(t)))
        e_b = end - starter
        e_e = end + starter

        end_l = len(t[e_b:e_e])
        #signs = rng.choice([-1,1], size=2*fft_bounds)
        #scale = signs*rng.exponential(scale=1, size=2*fft_bounds)
        #tau[(start - fft_bounds) : (start + fft_bounds)] += scale
        signs = rng.choice([-1,1], size=end_l)
        scale = signs*rng.exponential(scale=1, size=end_l)
        tau[(e_b) : e_e] += scale[::-1]

    return tau

def error(fit, solution, flag, t):
    if(flag == 0):
        err = np.abs(solution - fit)/np.abs(solution)
    else:
        err = np.sqrt(np.trapezoid(((fit - solution)**2), x=t)) / np.sqrt(np.trapezoid((solution**2), x=t))

    return err

    
