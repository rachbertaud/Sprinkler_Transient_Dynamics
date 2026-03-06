import matplotlib.pyplot as plt

def plot_analytical(full_t, full_y, y_fit, index):
    fig, ax = plt.subplots()

    # data as points
    ax.scatter(full_t, full_y, color='hotpink', s=10, label='data', zorder=2)

    # analytical fit as line
    ax.plot(full_t, y_fit, color='mediumpurple', linewidth=2, label='fit', zorder=3)

    # lime green circle around index
    ax.scatter(full_t[index], full_y[index], s=200, facecolors='none', 
               edgecolors='limegreen', linewidths=2.5, zorder=4, label='fit start')

        # zoom around fit region
    ax.set_xlim(full_t[index] - 1, full_t[-1] + 1)
    y_region = full_y[index:]
    ax.set_ylim(y_region.min() - 0.1, y_region.max() + 0.1)

    ax.legend()
    ax.set_xlabel('t')
    ax.set_ylabel('y')
    plt.tight_layout()
    plt.show()


def plot_franken(franken_t, franken_y, full_t, index, t_end):
    checker = False
    if(len(t_end) == 2048):
        checker = True
        n_neg = 2048
    else:
        n_neg = len(franken_t) - len(full_t[:index]) - len(t_end)
        

    t_zeros = franken_t[:n_neg]
    
    if checker:
        t_orig  = franken_t[n_neg:n_neg + index]
        t_tail  = franken_t[n_neg + index:]
    else:
        t_orig = franken_t[2048:4096]

    y_zeros = franken_y[:n_neg]
    
    if checker:
        y_orig  = franken_y[n_neg:n_neg + index]
        y_tail  = franken_y[n_neg + index:]
    else:
        y_orig = franken_y[2048:4096]

    plt.figure()
    plt.plot(t_zeros, y_zeros, color='cornflowerblue', linewidth=2, label='added zeros')
    plt.plot(t_orig,  y_orig,  color='hotpink',        linewidth=2, label='original data')
    if checker:
        plt.plot(t_tail,  y_tail,  color='mediumpurple',   linewidth=2, label='analytical tail')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y')
    plt.tight_layout()
    plt.show()


def plot_phi_gen(full_t, full_y, phi_gen):
    fig, ax = plt.subplots()

    ax.scatter(full_t, full_y, color='hotpink', s=10, label='data', zorder=2)
    ax.plot(full_t, phi_gen, color='mediumpurple', linewidth=2, label='phi_gen', zorder=3)

    ax.set_xlabel('t')
    ax.set_ylabel('y')
    ax.legend()
    plt.tight_layout()
    plt.show()
