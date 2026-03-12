import os
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_data(t_data, y_data, flag):
    # defines figure
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    # plots data
    line, = ax.plot(t_data, y_data, color='pink')
    vline = ax.axvline(x=t_data[0], color='purple', linestyle='--')

    # add slider
    ax_slider = plt.axes([0.2, 0.08, 0.6, 0.03])
    slider = Slider(ax_slider, 'Time of Choice', t_data.min(), t_data.max(), valinit=t_data[0], color="purple")

    def update(val):
        vline.set_xdata([slider.val, slider.val])
        fig.canvas.draw_idle()

    slider.on_changed(update)
    if(flag == 1):
        plt.title("Please use slider to pick time to start data fit. Exit plot once done.")
    if(flag == 2):
        plt.title("Please use slider to pick time to end data fit. Exit plot once done.")
    else:
        plt.title("Please use slider to pick time to end quiet region for noise characterization. Exit plot once done.")
    plt.show()

    return slider.val  # returns the final value when window is closed


def fit_segments(full_t, full_y, t_peaks, t_target, t_insert):
    peak_index = np.argmin(np.abs(t_peaks - t_target))
    insert_index = np.argmin(np.abs(full_t - t_insert))
    index = np.argmin(np.abs(full_t - t_peaks[peak_index]))
    
    t_seg = full_t[index:insert_index]
    y_seg = full_y[index:insert_index]
    
    return index, peak_index, insert_index, t_seg, y_seg 
