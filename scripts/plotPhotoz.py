#!/usr/bin/env python

"""
Jean coupon - 2015
script to plot results for project Stacked_X_ray
"""

import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages

from scipy import interpolate
#from scipy import linalg

from   astropy.io        import ascii,fits

import plot_utils
from scipy.optimize import curve_fit

from   astropy.cosmology import FlatLambdaCDM

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #

cosmo = FlatLambdaCDM(H0=72.0, Om0=0.258)

h = cosmo.H0.value/100.0

EPS    = np.finfo(np.float).eps
COLOR  = [ 'blue', 'green', 'red', 'orange', \
           'magenta', 'purple', 'lightblue', \
           'pink', 'Cyan', 'Brown', \
           'DarkRed', 'Indigo', 'DarkGreen', 'Salmon']

FIGX = 6.0
FIGY = 4.0

MARKER_SIZE = 2
FONT_SIZE = 16

np.random.seed( seed = 20091982)


BINS_TOMO     = np.array([0.2, 0.45, 0.55, 0.70, 0.80, 0.9, 1.0, 1.15, 1.35, 1.65, 2.0])
BINS_TOMO_CEN = 0.5 * (BINS_TOMO[1:] + BINS_TOMO[:-1])

STEP          = 1.e-3


# ----------------------------------------------------- #
# main
# ----------------------------------------------------- #

def main(args):
    function = getattr(sys.modules[__name__], args.option)(args, show=False)
    return


# ----------------------------------------------------- #
# Main functions
# ----------------------------------------------------- #



def plotPDF(args, show=False):

    fileInName=args.input.split(",")
    # PDF_keys = args.PDF_keys.split(",")

    if args.weight_key is not None:
        [zp, zs, weight], _ = getCols(fileInName[0], [args.zp_key, args.zs_key, args.weight_key], selection=args.select)
    else:
        [zp, zs], _ = getCols(fileInName[0], [args.zp_key, args.zs_key], selection=args.select)
        weight = None

    PDF, _ = getCols(fileInName[0], ['PDF'], selection=args.select)
    PDF = PDF[0]
    # PDF_bins = [float(s) for s in PDF_keys]


    with fits.open(fileInName[0]) as data:
        PDF_bins = data[2].data['BINS']


    # normalise PDF
    for i in range(len(PDF)):
    # for i in range(100):
        norm = int_trapz(PDF_bins, PDF[i], 0.0, 6.0)
        PDF[i] /= norm

        # to fix
        PDF[i,0] = 0.0

    if weight is not None:
        from scipy import stats

        [mode, _] = stats.mode(weight)
        weight /= mode

        indices=[]
        for i,w in enumerate(weight):
            for n in range(np.random.poisson(w)):
                indices.append(i)

        zp = zp[indices]
        zs = zs[indices]

        PDF = PDF[indices,:]


    #fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = 2

    #pp = PdfPages("test.pdf")
    pp = PdfPages(args.output)

    if True:
        mean  = np.zeros(len(BINS_TOMO_CEN))
        f_05  = np.zeros(len(BINS_TOMO_CEN))
        f_15  = np.zeros(len(BINS_TOMO_CEN))
        Ngals = np.zeros(len(BINS_TOMO_CEN))
        z     = np.arange(BINS_TOMO[0], BINS_TOMO[-1], STEP)
        for b in range(len(BINS_TOMO_CEN)):
            select = ((BINS_TOMO[b] < zp) & (zp < BINS_TOMO[b+1]))
            Ngals[b]  = len(zs[select])
            mean[b], f_05[b], f_15[b] = plot_PDF_dist(PDF_bins, PDF[select], zp[select], zs[select], BINS_TOMO_CEN[b], get_stats=True)


        plt.figure(figsize=(11, 10))
        plt.subplot(3,1,1); plt.ylabel('$\langle \Delta_z \\rangle$')
        plt.fill_between(z, -0.002+0.0*z, 0.002+0.0*z, facecolor='red', alpha = 0.5)
        plt.errorbar(BINS_TOMO_CEN, mean, yerr=mean/np.sqrt(Ngals), fmt='o')
        p = plt.Rectangle((0, 0), 1, 1, fc="r", alpha=0.5)
        plt.legend([p], ["Euclid requirements"], frameon=False, loc='upper left')


        plt.subplot(3,1,2); plt.ylabel('$f_{{0.05}}(\%)$'); plt.ylim(30.0 , 100.0)
        plt.errorbar(BINS_TOMO_CEN, 100.0*f_05, yerr=100.0*f_05/np.sqrt(Ngals), fmt='o')
        plt.fill_between(z, 68.0+0.0*z, 100.0+0.0*z, facecolor='red', alpha = 0.5)

        plt.subplot(3,1,3); plt.xlabel('$z_\mathrm{phot}}$');  plt.ylabel('$f_{{0.15}}(\%)$'); plt.ylim(70.0, 100.0)
        plt.errorbar(BINS_TOMO_CEN, 100.0*f_15, yerr=100.0*f_15/np.sqrt(Ngals), fmt='o')
        plt.fill_between(z, 90.0+0.0*z, 100.0+0.0*z, facecolor='red', alpha = 0.5)

        pp.savefig()


    if True:
        plt.figure(figsize=(11, 10))
        plt.subplot(2,1,1)
        for b in range(len(BINS_TOMO)-1):
            select = ((BINS_TOMO[b] < zp) & (zp < BINS_TOMO[b+1]))
            plot_nz(PDF_bins, PDF[select], zs[select], info=False, color=COLOR[b], xmax=2.5)

        plt.subplot(2,1,2)
        for b in range(len(BINS_TOMO)-1):
            select = ((BINS_TOMO[b] < zp) & (zp < BINS_TOMO[b+1]))
            plot_nz(PDF_bins, PDF[select], zs[select], info=False, true_nz=True, color=COLOR[b], xmax=2.5)


        pp.savefig()

        #plt.savefig("test.pdf")


    if True:

        for b in range(len(BINS_TOMO)-1):
            plt.figure(figsize=(11, 5))
            select = ((BINS_TOMO[b] < zp) & (zp < BINS_TOMO[b+1]))
            plt.subplot(1,2,1)
            plt.title('${0:3.3f} < z_\mathrm{{phot}} < {1:3.3f}$'.format(BINS_TOMO[b], BINS_TOMO[b+1]))
            plot_nz(PDF_bins, PDF[select], zs[select])
            plt.subplot(1,2,2)
            plot_PDF_dist(PDF_bins, PDF[select], zp[select], zs[select], BINS_TOMO_CEN[b])

            pp.savefig()





    pp.close()


    return


def plotSimuSigma(args, show=False):

    fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = 2
    plot_utils.axes(ax, "$\mathrm{Flux\,error}_\mathrm{DC2} \\times $" , "$\sigma/(1+z)$", [1.70, 0.01], [0.03, 0.08], title=args.title)

    data  = ascii.read("{0:s}".format(args.input), format="commented_header", header_start=-1)

    #N = [float(n) for n in data['N']]
    err = 1.0/np.sqrt(12000.0)*data['sigma']

    plot_utils.markers(ax, data['factor'], data['sigma'], err, size, COLOR[0], 'Simulation')

    popt, pcov = curve_fit(power_law, data['factor'], data['sigma'], sigma=err)

    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1])

    ax.plot(x, power_law(x, popt[0], popt[1]), color= COLOR[0], lw=2)

    plot_utils.markers(ax, [1.0], [0.057], err[0], size, COLOR[1], 'DC2')

    ax.plot(ax.get_xlim(), [0.05, 0.05], color=COLOR[2], label='Requirement', lw=size)

    #return
    # ----------------------------------------------------------------- #
    # save figure
    # ----------------------------------------------------------------- #

    plt.legend(loc='upper right', frameon=False, numpoints=1)


    fig.set_tight_layout(True)
    fig.savefig(args.output)

    if show:
        plt.show()

    return

def plotSimuEta(args, show=False):

    fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = 2
    plot_utils.axes(ax, "$\mathrm{Flux\,error}_\mathrm{DC2} \\times $" , "$\eta$", [1.70, 0.01], [5.0, 20.0], title=args.title)

    data  = ascii.read("{0:s}".format(args.input), format="commented_header", header_start=-1)

    #N = [float(n) for n in data['N']]
    err = 1.0/np.sqrt(12000.0)*data['eta']

    plot_utils.markers(ax, data['factor'], data['eta'], err, size, COLOR[0], 'Simulation')

    popt, pcov = curve_fit(power_law, data['factor'], data['eta'], sigma=err)

    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1])

    ax.plot(x, power_law(x, popt[0], popt[1]), color= COLOR[0], lw=2)

    plot_utils.markers(ax, [1.0], [15.64], err[0], size, COLOR[1], 'DC2')

    ax.plot(ax.get_xlim(), [10.0, 10.0], color=COLOR[2], label='Requirement', lw=size)

    #return
    # ----------------------------------------------------------------- #
    # save figure
    # ----------------------------------------------------------------- #

    plt.legend(loc='upper right', frameon=False, numpoints=1)


    fig.set_tight_layout(True)
    fig.savefig(args.output)

    if show:
        plt.show()

    return



def power_law(x, a, b):
    return a * x + b


def plotScatter(args, show=False):

    zmin = args.zmin
    zmax = args.zmax

    fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = FONT_SIZE

    plot_utils.axes(ax, "$z_\mathrm{spectro}$" , "$z_\mathrm{photo}$", [zmin, zmax], [zmin, zmax], title=args.title)

    # ----------------------------------------------------------------- #
    # scatter plot
    # ----------------------------------------------------------------- #

    fileInName=args.input.split(",")

    if args.weight_key is not None:
        [zp, zs, weight], _ = getCols(fileInName[0], [args.zp_key, args.zs_key, args.weight_key], selection=args.select)
    else:
        [zp, zs], _ = getCols(fileInName[0], [args.zp_key, args.zs_key], selection=args.select)
        weight = None

    if weight is not None:
        from scipy import stats

        # [mode, _] = stats.mode(weight)
        # weight /= mode

        indices=[]
        for i,w in enumerate(weight):
            for n in range(np.random.poisson(w)):
                indices.append(i)

        zp = zp[indices]
        zs = zs[indices]

    if args.density:
        hh, locx, locy = np.histogram2d(zs, zp, range=[[zmin, zmax],[zmin, zmax]], bins=[150, 150])
        # hh, locx, locy = np.histogram2d(zs, zp, range=[[zmin, zmax],[zmin, zmax]], bins=[80, 80])
        hh[hh < EPS] = np.nan
        ax.imshow(np.log(hh.T), origin="lower", cmap='Blues', extent=np.array([[zmin, zmax],[zmin, zmax]]).flatten(), aspect='auto', interpolation="nearest")
    else:
        plot_utils.markers(ax, zs, zp, None, 0.3, COLOR[0], "")

    # ----------------------------------------------------------------- #
    # stats
    # ----------------------------------------------------------------- #

    sigma, eta, bias, eta2sig = stats_zpzs(zp, zs, [zmin, zmax])

    stats_string  = '$N_\mathrm{{gals}}$ = {0}'.format(len(zp))
    stats_string += '\n$\sigma = {0:5.3f} \\times (1+z)$'.format(sigma[0])
    stats_string += '\n$\eta = {0:5.2f}\%$'.format(eta[0])
    stats_string += '\nbias$ = {0:5.3f} \\times (1+z)$'.format(bias[0])
    # stats_string += '\n$\eta_{{2\sigma}}/(1+z) = {0:5.2f}\%$'.format(eta2sig[0])

    ax.text(0.05*(zmax-zmin)+zmin, 0.75*(zmax-zmin)+zmin, stats_string, size=FONT_SIZE)

    # ----------------------------------------------------------------- #
    # red lines
    # ----------------------------------------------------------------- #

    x = np.arange(zmin, zmax, 0.01)
    ax.plot(x, x, 'r-', lw=2)
    ax.plot(x, x + 0.15*(1+x), 'r:', lw=2)
    ax.plot(x, x - 0.15*(1+x), 'r:', lw=2)

    # ----------------------------------------------------------------- #
    # save figure
    # ----------------------------------------------------------------- #

    fig.set_tight_layout(True)
    fig.savefig(args.output)

    if show:
        plt.show()

    return



    #print len(zp)

    return

# ----------------------------------------------------- #
# Utils
# ----------------------------------------------------- #


def stats_zpzs(zp, zs, bins):

    sigma   = np.zeros((len(bins)-1))
    eta     = np.zeros((len(bins)-1))
    eta2sig = np.zeros((len(bins)-1))
    bias    = np.zeros((len(bins)-1))

    for i in range(0, len(bins)-1):
        select          = ((bins[i] < zp) & (zp < bins[i+1])).nonzero()
        dist            = (zp[select]-zs[select])/(1+zs[select])
        n_outliers      = (abs(dist) > 0.15).nonzero()

        if(len(dist) > 0):
            sigma[i]          = 1.48*np.median(np.abs(dist))
            eta[i]            = float(len(dist[n_outliers]))/float(len(dist)) * 100.0
            bias[i]           = np.median(dist)
            n_outliers_2sig   = (abs(dist) > 2.0*sigma[i]).nonzero()
            eta2sig[i]        = float(len(dist[n_outliers_2sig]))/float(len(dist)) * 100.0

    return sigma, eta, bias, eta2sig


def getCols(fileInName, keys, selection=None, array=False):
    '''
    Return column arrays from fileInName (fits file)

    INPUT
    fileInName: fits file name
    keys: column names list or expression (e.g. "col1 - col2")
    selection: selection (string)
    array: if set returns 2D array (N,K)

    OUTPUT
    list of column arrays
    (optional) selection array

    WARNING
    Numeric column names cannot be used
    '''

    # first open the fits file and read the data
    fileIn = fits.open(fileInName)
    data = fileIn[1].data
    cols = fileIn[1].columns
    fileIn.close()

    if selection is not None:
        cmd = ""
        for s in re.split('(\W+)', selection):
            if s in cols.names:
                cmd += "data['"+s+"']"
            else:
                cmd += s

        select = eval(cmd)
        str_select = "[select]"
    else:
        str_select = ""

    result=[]
    # replace key value with command
    for k in keys:
        cmd = ""
        for s in re.split('(\W+)', k):
            if s in cols.names:
                cmd += "data"+str_select+"['"+s+"']"
            else:
                cmd += s

        # print cmd
        result.append(eval(cmd))

    if array:
        K = len(result)
        N = len(result[0])
        result_tmp = np.zeros((N,K))
        for d in range(K):
            result_tmp[:,d] = result[d]

        result = result_tmp

    # return result
    if selection is not None:
        return result, select
    else:
        return result, None


def int_trapz(x, y, a, b):
    # return the integrated value between a and b
    # use the trapeze rule
    int_range = np.linspace(a, b, num=1000)
    f         = np.interp(int_range, x, y)
    return np.sum(np.diff(int_range) * (f[:-1]+f[1:])/2.0)


def plot_PDF_dist(PDF_bins, PDF, zphot, zspec, z0, get_stats=False):

    # Compute the PDF(z-z_true) sum
    PDF_dist, PDF_dist_bins = get_PDF_dist(PDF_bins, PDF, -4.0, +4.0, zspec)
    #hist, hist_bins = np.histogram(zphot-zspec, bins=50, density=True)
    #PDF_dist = hist
    #PDF_dist_bins = 0.5*(hist_bins[1:]+hist_bins[:-1])


    # correct for the bias as well ??

    # mode, mean and standard deviation of the distribution
    mode    = PDF_dist_bins[np.argmax(PDF_dist)]
    mean    = int_trapz(PDF_dist_bins, PDF_dist_bins*PDF_dist, -4.0, 4.0)/(1.0+z0)
    std_dev = np.sqrt(int_trapz(PDF_dist_bins, np.square(PDF_dist_bins-mean) * PDF_dist, -4.0, 4.0))/(1.0+z0)

    # fraction of probability within 0.05(1+z) and 0.15(1+z)
    z_68 = 0.05*(1.0+z0)
    z_99 = 0.15*(1.0+z0)

    f_05 = int_trapz(PDF_dist_bins-mode, PDF_dist, -z_68, +z_68)
    f_15 = int_trapz(PDF_dist_bins-mode, PDF_dist, -z_99, +z_99)

    if get_stats:
        return mean, f_05, f_15

    ax = plt.plot(PDF_dist_bins, PDF_dist, color='b', lw=3)
    plt.fill_between(PDF_dist_bins, 0.0, PDF_dist, color='b', alpha=0.5, label='$\Sigma$',\
                     where= (-z_68+mode < PDF_dist_bins) & (PDF_dist_bins < +z_68+mode))
    plt.fill_between(PDF_dist_bins, 0.0, PDF_dist, color='b', alpha=0.5, label='$n(z_{{spec}}$',\
                     where= (-z_99+mode < PDF_dist_bins) & (PDF_dist_bins < +z_99+mode))

    xlim = plt.xlim(-2.0, 2.0)
    ylim = plt.ylim(0.0,); yspan = ylim[1]-ylim[0]

    plt.xlabel('$\Delta z$')
    plt.ylabel('PDF$(z-z_{\mathrm{spec}})$')

    plt.axvline(x=mean, color='b', lw=3)
    plt.axvline(x=0.0, color='r')

    plt.text(xlim[0], 0.90*yspan + ylim[0], ' $\mathrm{{[in\,(1+z)\,unit]}}$', fontsize=15)
    plt.text(xlim[0], 0.83*yspan + ylim[0], ' $f_{{0.05}} = {0:3.2f}\%$'.format(f_05*100), fontsize=15)
    plt.text(xlim[0], 0.76*yspan + ylim[0], ' $f_{{0.15}} = {0:3.2f}\%$'.format(f_15*100), fontsize=15)
    plt.text(xlim[0], 0.69*yspan + ylim[0], ' $\langle \Delta_z \\rangle = {0:3.3f}$'.format(mean), fontsize=15)
    #plt.text(xlim[0], 0.62*yspan + ylim[0], ' $\langle \Delta_z^2 \\rangle = {0:3.3f}$'.format(std_dev), fontsize=15)



def get_PDF_dist(x, PDF, a, b, x0=None):
    # returns the sum of normalised PDF(x-x0)
    # PDF should have linearly spaced values

        Nobjects  = len(PDF)

        int_range = np.arange(a, b, STEP)
        result    = np.zeros(len(int_range))

        # loop over objects
        if x0 is None:
            for i in range(Nobjects):
                if np.sum(PDF[i]) > 0.0:
                    result += np.interp(int_range, x, PDF[i]/np.sum(PDF[i])/(x[1] - x[0]))
        else:
            for i in range(Nobjects):
                if np.sum(PDF[i]) > 0.0:
                    result += np.interp(int_range, x-x0[i], PDF[i]/np.sum(PDF[i])/(x[1] - x[0]))

        # return normalised sum and bins
        return result/Nobjects, int_range

def plot_nz(PDF_bins, PDF, zspec, info=True, true_nz=False, color='b', xmax=4.0):

    # plot options
    plt.xlabel('$z$')
    plt.ylabel('$n(z)$')
    xlim = plt.xlim(0.0, xmax)

    # Compute the PDF(z) sum
    if true_nz :
        hist, hist_bins = np.histogram(zspec, bins=50, density=True)
        plt.fill_between(0.5*(hist_bins[1:]+hist_bins[:-1]), 0.0, hist, color=color, alpha=0.5)
        plt.plot(0.5*(hist_bins[1:]+hist_bins[:-1]), hist, color=color, lw=2)
        a = plt.Rectangle((0, 0), 1, 1, fc="Silver", ec="Silver", lw = 2, alpha=0.5)
        plt.legend([a], ['$n(z_\mathrm{zspec})$'], frameon=False, loc='center right')
    # Compute the true n(z)
    else:
        PDF_nz, PDF_nz_bins = get_PDF_dist(PDF_bins, PDF, +0.0, +4.0)
        plt.fill_between(PDF_nz_bins, 0.0, PDF_nz, color=color, alpha=0.5)
        plt.plot(PDF_nz_bins, PDF_nz, color=color, lw=2)
        a = plt.Rectangle((0, 0), 1, 1, fc="Silver", ec="Silver", lw = 2, alpha=0.5)
        plt.legend([a], ['$\Sigma \mathrm{PDF}$'], frameon=False, loc='center right')

    # plot true n(z) plus info
    if info:
        hist, hist_bins = np.histogram(zspec, bins=50, density=True)
        plt.fill_between(0.5*(hist_bins[1:]+hist_bins[:-1]), 0.0, hist, color='r', alpha=0.5)
        plt.plot(0.5*(hist_bins[1:]+hist_bins[:-1]), hist, color='r', lw=2)

        Ngals = len(PDF)
        mean  = int_trapz(PDF_nz_bins, PDF_nz_bins*PDF_nz, -4.0, 4.0)

        a = plt.Rectangle((0, 0), 1, 1, fc="b", ec="b", lw = 2, alpha=0.5)
        b = plt.Rectangle((0, 0), 1, 1, fc="r", ec="r", lw = 2, alpha=0.5)
        plt.legend([a, b], ['$\Sigma \mathrm{PDF}$','$n(z_\mathrm{zspec})$'], frameon=False, loc='center right')

        ylim = plt.ylim(0.0, ); yspan = ylim[1]-ylim[0]
        plt.text(2.0, 0.90*yspan+ylim[0], '$n_\mathrm{{gals}} = {0}$'.format(Ngals), fontsize=15)
        plt.text(2.0, 0.83*yspan+ylim[0], '$\langle z_\mathrm{{PDF}}\\rangle = {0:3.3f}$'.format(mean), fontsize=15)
        plt.text(2.0, 0.76*yspan+ylim[0], '$\langle z_\mathrm{{spec}}\\rangle = {0:3.3f}$'.format(np.mean(zspec)), fontsize=15)




# ----------------------------------------------------- #
# Main
# ----------------------------------------------------- #

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('option',                                   help="Which quantity to plot")
    parser.add_argument('-i', '--input',              default=None, help='input file')
    parser.add_argument('-o', '--output',             default=None, help='output file')
    parser.add_argument('-t', '--title',  type=str,   default="",   help='Title')
    parser.add_argument('-s', '--select', type=str,   default=None,   help='Selection')
    parser.add_argument('-zp_key', type=str,   default="Z_BEST",    help='Photo-z key')
    parser.add_argument('-zs_key', type=str,   default="ZSPEC",     help='Spec-z key')
    parser.add_argument('-weight_key', type=str,   default=None,     help='Weight key')
    parser.add_argument('-density', action='store_true', help='Display point density')
    parser.add_argument('-zmax',   type=float, default=4.0,  help='Upper redshift limit')
    parser.add_argument('-zmin',   type=float, default=0.0,  help='Lower redshift limit')

    args = parser.parse_args()

    main(args)
