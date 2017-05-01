from params import *

def timer(start, end, label):
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:05.2f}".format(int(minutes), seconds), '\t', label)


def plotGrid():
    for i in range(int(N_x * n_subdiv)):
        if i % n_subdiv == 0:
            plt.axvline(cell_size * i / n_subdiv, color='r')
        else:
            plt.axvline(subcell_size * i, color='b', linestyle='dotted')
    for j in range(int(N_y * n_subdiv)):
        if j % n_subdiv == 0:
            plt.axhline(cell_size * j / n_subdiv, color='r')
        else:
            plt.axhline(subcell_size * j, color='b', linestyle='dotted')
    plt.axis((0, X_max, 0, Y_max))
    plt.show()


def plotDistribution(Array, map_bounds, label, interpolation, line):
    cmap = plt.cm.jet  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]  # extract all colors from the .jet map
    cmaplist[0] = (.5, .5, .5, 1.0)  # force the first color entry to be grey
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)  # create the new map
    norm = mpl.colors.BoundaryNorm(map_bounds, cmap.N)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(Array.T, origin='lower', cmap=cmap, norm=norm, interpolation=interpolation)
    fig.colorbar(im, ax=ax)

    if label == 'Subcells':
        xtick_locs = [i * 20 for i in range(int(len(Array[:, 0]) / 20) + 1)]
        ytick_locs = [i * 20 for i in range(int(len(Array[0, :]) / 20) + 1)]
        xtick_lbls = [i * 10 for i in range(int(N_x * n_subdiv / 20) + 1)]
        ytick_lbls = [i * 10 for i in range(int(N_x * n_subdiv / 20) + 1)]
    else:
        xtick_locs = [i * 10 for i in range(int(len(Array[:, 0]) / 10) + 1)]
        ytick_locs = [i * 10 for i in range(int(len(Array[0, :]) / 10) + 1)]
        xtick_lbls = [i * 10 for i in range(int(N_x / 10) + 1)]
        ytick_lbls = [i * 10 for i in range(int(N_x / 10) + 1)]

    plt.xticks(xtick_locs, xtick_lbls)
    plt.yticks(ytick_locs, ytick_lbls)
    plt.title(label)
    if line:
        plt.plot([i for i in range(N_x)], [i for i in range(N_y)], color='r')
    plt.show()
    fig.clf()
    plt.close()

    del ax, im, fig, cmap, norm
    gc.collect()

def plot_T_function(Array):
    def f(x,a,b):
        return a*x+b

    T = np.zeros((N_x))
    for i in range(N_x):
        T[i] = Array[i, i]
    x = np.array([i for i in range(N_x)])
    popt, pcov = opt.curve_fit(f, x, T)

    plt.plot(T, label = 'actual')
    plt.plot(x, f(x, *popt), 'r--', label='fitted')
    plt.xlabel('cells')
    plt.ylabel('T')
    plt.legend(loc=0)
    plt.show()

def plot_N_function(Array):
    def f(x,a,b):
        return a*x+b

    N = np.zeros((N_x))
    for i in range(N_x):
        N[i] = Array[i, i]

    x = np.array([i for i in range(N_x)])
    popt, pcov = opt.curve_fit(f, x, N)

    plt.plot(N, label = 'actual')
    plt.plot(x, f(x, *popt), 'r--', label='fitted')
    plt.xlabel('cells')
    plt.ylabel('N')
    plt.legend(loc=0)
    plt.show()

def plotVelDistribusion(beta):
    bin_num = 50

    def gauss(u, beta):
        return beta / sqrt(pi) * np.exp(-beta ** 2 * u ** 2)

    def speed_pdf(c, beta):
        return 4 / sqrt(pi) * beta ** 3 * c ** 2 * np.exp(-beta ** 2 * c ** 2)

    def fitGaussian(field):
        x = np.arange(min(field), max(field), 10)
        plt.hist(field, bins=bin_num, normed=1, alpha=0.4)
        plt.plot(x, gauss(x, beta), 'r--', linewidth=3, label=r'$\frac{\beta}{\sqrt{\pi}}\exp(-\beta^2u^2)$')
        plt.legend(loc=0)

    # U_x pdf
    # field = Vel[:,0]
    # fitGaussian(field)
    # plt.xlabel(r'$u^{\prime}$', fontsize='xx-large')
    # plt.ylabel(r'$f_{u^{\prime}}$', fontsize='xx-large')
    # plt.show()
    # # U_y pdf
    # field = Vel[:,1]
    # fitGaussian(field)
    # plt.xlabel(r'$v^{\prime}$', fontsize='xx-large')
    # plt.ylabel(r'$f_{v^{\prime}}$', fontsize='xx-large')
    # plt.show()
    # # U_z
    # field = Vel[:,2]
    # fitGaussian(field)
    # plt.xlabel(r'$w^{\prime}$', fontsize='xx-large')
    # plt.ylabel(r'$f_{w^{\prime}}$', fontsize='xx-large')
    # plt.show()
    # speed pdf
    field = np.sqrt(Vel[:, 0] ** 2 + Vel[:, 1] ** 2 + Vel[:, 2] ** 2)
    x = np.linspace(min(field), max(field), 100)
    plt.hist(field, bins=bin_num, normed=1, alpha=0.4)
    plt.plot(x, speed_pdf(x, beta), 'r--', linewidth=3,
             label=r'$\frac{4}{\sqrt{\pi}}\beta^3c^2exp(-\beta^2c^{\prime 2})$')
    plt.xlabel(r'$c^{\prime}$', fontsize='xx-large')
    plt.ylabel(r'$f_{c^{\prime}}$', fontsize='xx-large')
    plt.title(r'Histogram with ' + str(len(Vel[:, 0])) + " samples")
    plt.legend(loc=0)
    plt.show()
