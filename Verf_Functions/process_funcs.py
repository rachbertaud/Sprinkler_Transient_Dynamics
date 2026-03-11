import numpy as np

def combine_data(full_t, full_y, index, phi_an, proc_data_switch, N_f):
    half = N_f//2
    
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
        t_neg = np.arange(-half * dt, 0, dt)
        y_neg = np.zeros(len(t_neg))
    else:
        t_neg = np.arange(-half * dt, 0, dt)
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
