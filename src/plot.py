import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
# import time
# import progressbar as pgbar
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc, rcParams
import matplotlib.colors as colors

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=10)
plt.rcParams['savefig.dpi'] = 200

def plot_dens_profs():
    inputadr = "../lh1.out"

    # Input file:
    if not os.path.isfile(inputadr):
        raise IOError()
    print("\nReading data...")
    # inputf = open(inputadr, "r")
    # data = ascii.read(inputf)
    data_tbl = np.loadtxt(inputadr).T
    print("Done!\n")

    colnames = ['istep',
                'R',
                'ix',
                'Rsh',
                'GammaSh',
                'rhoSh',
                'rho(ix)'
                ]

    # Variables to plot:
    sx = 'R'
    sy = 'rho(ix)'
    xcol = data_tbl[colnames.index(sx)]
    ycol = data_tbl[colnames.index(sy)]

    i = 0
    finished = False
    startPos = 0
    f, axes = plt.subplots(ncols=1,nrows=2, sharex = "all")
    cmap = cm.Blues
    norm = Normalize(vmin=0, vmax=len(set(data_tbl[colnames.index('istep')])))

    while not(finished):
        istep = data_tbl[colnames.index('istep')][startPos]
        condi = (data_tbl[colnames.index('istep')] == istep)
        x = xcol[condi]
        y = ycol[condi]

        ymax = np.max(ycol[~np.isnan(ycol)])
        ymin = np.min(ycol[~np.isnan(ycol)])
        xmax = np.max(xcol[~np.isnan(xcol)])
        xmin = np.min(xcol[~np.isnan(xcol)])
        startPos += len(x)

        ax = axes[0]
        ax.set_ylim(ymin,ymax)
        ax.set_ylabel(r"$\rho$")
        ax.plot(x,y,color=cmap(norm(int(istep))))
        ax.set_xscale("log")
        ax.set_yscale("log")

        ax = axes[1]
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_ylabel(r"$\rho$")
        ax.plot(x,y,color=cmap(norm(int(istep))))
        ax.set_xscale("log")
        ax.set_yscale("log")

        if startPos == len(data_tbl[0]):
            finished = True
        i += 1

    plt.xlabel("Radius")
    plt.show()

def main():
    plot_dens_profs()

if __name__ == '__main__':
    main()