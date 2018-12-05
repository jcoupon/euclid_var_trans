#!/usr/bin/env python

"""
Jean coupon - 2017
scripts to compute photo-z PDF from NNPZ
"""

import os, sys
import numpy as np
import re

from scipy import spatial
from scipy import interpolate
from scipy.stats import gaussian_kde
import collections
from astropy.io import ascii,fits

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #

EPS    = np.finfo(np.float).eps
COLOR  = [ 'blue', 'green', 'red', 'orange', \
           'magenta', 'purple', 'lightblue', \
           'pink', 'Cyan', 'Brown', \
           'DarkRed', 'Indigo', 'DarkGreen', 'Salmon']

FIGX = 6.0
FIGY = 4.0

MARKER_SIZE = 2
FONT_SIZE = 16

# np.random.seed( seed = 20091983)

# ----------------------------------------------------- #
# main
# ----------------------------------------------------- #

def main(args):
    function = getattr(sys.modules[__name__], args.option)(args)
    return

# ----------------------------------------------------- #
# Main functions
# ----------------------------------------------------- #

def getWeight(args):
    '''
    Compute the weight for a calibration file

    INPUT
    - fileInNames: fits file names: calibration_file,reference_file
    - keys: column names describing dimensions, example; "MAG_I","MAG_G-MAG_R","MAG_R-MAG_I"
    - selection: selection string for both files

    NOTES:
    - both the photoz_file and reference_file must contain
    the same columns as defined by "keys"

    OUTPUT
    - calibration weights
    '''

    verbose = True
    Nnei = 100

    fileInName = args.input.split(",")
    selection = args.select.split(",")

    if verbose: sys.stderr.write("Reading input files...")

    # calib_file
    ref, select = getCols(fileInName[0], args.keys.split(","), selection=selection[0], array=True)
    (Nref, Ndim) = ref.shape
    ref[np.logical_not(np.isfinite(ref))] = 0.0 # set NaNs and inf to 0


    # sample file
    sample, select = getCols(fileInName[1], args.keys.split(","), selection=selection[1], array=True)
    (Nsample, Ndim) = sample.shape
    sample[np.logical_not(np.isfinite(sample))] = 0.0 # set NaNs and inf to 0

    if verbose: sys.stderr.write("done\n")

    if verbose: sys.stderr.write("Building trees...")
    refTree    = spatial.KDTree(ref)
    sampleTree = spatial.KDTree(sample)
    if verbose: sys.stderr.write("done\n")

    # return

    w = np.zeros(Nref)
    for i in range(Nref):
    #for i in range(20):

        (d, _) = refTree.query(ref[i,:], Nnei)
        indices = sampleTree.query_ball_point(ref[i,:], d[-1])

        w[i] = float(len(indices))

        if verbose:
            if (i+1)%100 == 0:
                sys.stderr.write("\r" + "Weights: computed {0:d} objects {1:5.2f}%".format(i+1, 100.0*float(i+1) / float(Nref)))
                sys.stderr.flush()
    if verbose: sys.stderr.write("\r" + "Weights: computed {0:d} objects {1:5.2f}%\n".format(i+1, 100.0))

    w *= float(Nref)/float(Nsample) / float(Nnei)

    if verbose: sys.stderr.write("Writing output file...")

    cols = []
    cols.append(fits.Column(name='weight', format='E', array=w))

    if args.merge_with_input:
        fileIn = fits.open(fileInName[0])
        hdu_0  = fileIn[0]
        hdu_1  = fits.BinTableHDU.from_columns(fileIn[1].columns + fits.ColDefs(cols))
        if len(fileIn) > 2:
            tbhdu  = fits.HDUList([hdu_0, hdu_1, fileIn[2]])
        else:
            tbhdu  = fits.HDUList([hdu_0, hdu_1])

    else:
        hdu_0 = fits.PrimaryHDU()
        hdu_1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        tbhdu  = fits.HDUList([hdu_0, hdu_1])
    tbhdu.writeto(args.output, overwrite=True)

    if args.merge_with_input:
        fileIn.close()

    if verbose: sys.stderr.write("done\n")



    return





def getPz(args):
    '''
    Compute the P(z) from nearest neighbors P(z)'s
    given a reference (training) file.

    INPUT
    - fileInNames: fits file names: inputFile^trainingFile
    - keys: column names describing dimensions, example; "MAG_I,MAG_G-MAG_R,MAG_R-MAG_I"
    - selection: selection strings inputSelect^trainingSelect

    NOTES:
    - both the photoz_file and reference_file must contain
    the same column names as defined by "keys"

    OUTPUT
    - P(z)'s of input file objects + input columns if merge_with_input is set
    '''

    """
    options
    """

    # for large dataset
    sys.setrecursionlimit(100000)

    from sklearn.neighbors import KernelDensity

    verbose = True

    fileInName=args.input.split("^")

    if args.select is not None:
        fileInSelect = args.select.split("^")
    else:
        fileInSelect = [None for f in fileInName]

    if verbose: sys.stderr.write("Reading input files...")

    """
    input file
    """
    sample, sampleSelect = getCols(fileInName[0], args.keys.split(","), selection=fileInSelect[0], array=True)

    if args.keys_err is not None:
        sample_err, _ = getCols(fileInName[0], args.keys_err.split(","), selection=fileInSelect[0], array=True)
    (Nsample, Ndim) = sample.shape
    sample[np.logical_not(np.isfinite(sample))] = 0.0 # set NaNs and inf to 0

    if Nsample == 0:
        raise ValueError("The input file is empty, exiting...")

    """
    training file
    """
    ref, select = getCols(fileInName[1], args.keys_ref.split(","), selection=fileInSelect[1], array=True)
    if args.keys_err is not None:
        ref_err, select = getCols(fileInName[1], args.keys_err.split(","), selection=fileInSelect[1], array=True)
    (Nref, Ndim) = ref.shape
    ref[np.logical_not(np.isfinite(ref))] = 0.0  # set NaNs and inf to 0


    # for reference: histo or sum of PDFs
    if args.PDF_histo:
        # if histogram, set bins
        bins = np.linspace(0.0, 6.0, num=601, endpoint=True)
        PDF_ref_bins = 0.5*(bins[1:]+bins[:-1])
        # get reference redshifts and weights
        # set with -keys_histo redshift,weight
        keys_histo = args.keys_histo.split(",")
        if len(keys_histo) > 1:
            #[z_ref, weight, source], select = getCols(fileInName[1], keys_histo, selection=fileInSelect[1])
            [z_ref, weight], select = getCols(fileInName[1], keys_histo, selection=fileInSelect[1])
        else:
            [z_ref], select = getCols(fileInName[1], keys_histo, selection=fileInSelect[1])
            weight = np.ones(len(z_ref))
    else:
        # if PDF, recover PDF from reference
        PDF_ref, PDF_ref_bins = getPDF(fileInName[1], normalise=args.no_norm, select=select, PDF_key="PDF_L15")
        weight = np.ones(len(PDF_ref))

    if verbose: sys.stderr.write("done\n")

    """
    Build tree of colors for reference
    """
    if verbose: sys.stderr.write("Building reference tree...")
    tree = spatial.KDTree(ref)
    if verbose: sys.stderr.write("done\n")

    # does not work:
    # import pickle
    # pickle.dump(tree, open("ref.pickle", 'w+'))
    # tree = pickle.load(open("ref.pickle", 'rb'))


    # return

    """
    main loop
    """
    # tests    Nsample = 1000

    Nnei = 50
    kde = False # see https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/

    # output arrays
    pofz = np.zeros((Nsample, len(PDF_ref_bins)), dtype=np.float32)
    if args.PDF_histo:

        zTrain = np.zeros((Nsample, Nnei), dtype=np.float32)
        if len(keys_histo) > 1:
            wTrain = np.zeros((Nsample, Nnei), dtype=np.float32)
            # sTrain = np.zeros((Nsample, Nnei), dtype=np.float32)

    est = collections.OrderedDict()
    for name in ['zmean','zmode','zmedian','z_std', 'zl95', 'zl68', 'zh68', 'zh95', 'z_mc', 'zconf']:
        est[name] = np.zeros(Nsample) -99.0

    zmin = PDF_ref_bins[0]
    zmax = PDF_ref_bins[-1]

    for i in range(Nsample):
    # for i in range(1):

        # photometry failure
        if abs(np.sum(sample[i,:])) < EPS:
            continue

        # find nearest neighbors
        (d, indices) = tree.query(sample[i,:], Nnei)

        # associated errors
        derr = np.ones(Nnei)
        if args.keys_err is not None:

            for n, j in enumerate(indices):
                err = np.sum(sample_err[i,:]+ref_err[j, :])
                if err > EPS:
                    derr[n] = err

        # weight: 1/distance * 1/err * input_weight
        w = np.ones(len(indices))
        for n, j in enumerate(indices):
            if d[n] > EPS:
                w[n] = 1.0/d[n] * 1.0/derr[n] * weight[j]

        if args.PDF_histo:


            zTrain[i, :] = z_ref[indices]
            if len(keys_histo) > 1:
                wTrain[i, :] = w
                # sTrain[i, :] = source[indices]

            if kde:

                z_weighted = weightedSample(z_ref[indices], w)

                std = np.std(z_weighted)

                if (std < 1.e-4) | (len(z_weighted) < 2):
                    pofz[i, :], _ = np.histogram(z_ref[indices], bins=bins, density=True, weights=w)
                else:
                    density    = gaussian_kde(z_weighted, bw_method=0.03/std)
                    pofz[i, :] = density.pdf(PDF_ref_bins)

                    # density = KernelDensity(kernel='gaussian', bandwidth=0.02).fit(z_weighted[:, np.newaxis])
                    # pofz[i, :] = density.score_samples(PDF_ref_bins[:, np.newaxis])


            else:
                pofz[i, :], _ = np.histogram(z_ref[indices], bins=bins, density=True, weights=w)
        else:
            for n, j in enumerate(indices):
                if (d[n] > EPS) & (np.sum(PDF_ref[j,:]) > EPS):
                    pofz[i, :] += PDF_ref[j,:] * w[n]

        # PDF_sample_inter = np.interp(PDF_ref_bins, PDF_sample_bins, PDF_sample[i,:])
        # pofz[i,:] *= PDF_sample_inter

        # normalise PDF
        norm = int_trapz(PDF_ref_bins, pofz[i,:],PDF_ref_bins[0], PDF_ref_bins[-1] )
        if norm > EPS:
            pofz[i,:] /= norm

            est['zmean'][i] = int_trapz(PDF_ref_bins, pofz[i,:]*PDF_ref_bins, zmin, zmax)
            est['zmode'][i] = max_pos_PDF(PDF_ref_bins, pofz[i,:])
            est['zmedian'][i] = medianfromDist(PDF_ref_bins, pofz[i,:])
            est['z_std'][i] = np.sqrt(int_trapz(PDF_ref_bins, pofz[i,:]*pow(PDF_ref_bins-est['zmean'][i], 2.0), zmin, zmax))
            est['zl95'][i] = sampleFromDist(PDF_ref_bins, pofz[i,:], q=0.05/2.0)
            est['zl68'][i] = sampleFromDist(PDF_ref_bins, pofz[i,:], q=0.32/2.0)
            est['zh68'][i] = sampleFromDist(PDF_ref_bins, pofz[i,:], q=1.0-0.32/2.0)
            est['zh95'][i] = sampleFromDist(PDF_ref_bins, pofz[i,:], q=1.0-0.05/2.0)
            est['z_mc'][i] = sampleFromDist(PDF_ref_bins, pofz[i,:])
            est['zconf'][i] = int_trapz(PDF_ref_bins, pofz[i,:], est['zmedian'][i]-0.03*(1.0+est['zmedian'][i]), est['zmedian'][i]+0.03*(1.0+est['zmedian'][i]))

        # for name in est.keys():
            # print "{0:s}:{1}".format(name,est[name][i])


            # test
            #z_weighted = weightedSample(PDF_ref_bins[pofz[i,:]>EPS], pofz[i,:][pofz[i,:]>EPS])
            #density = gaussian_kde(z_weighted)
            #pofz_KDE = density.pdf(PDF_ref_bins)
            #z_median[i] = medianfromDist(PDF_ref_bins, pofz_KDE)

        if verbose:
            if (i+1)%1000 == 0:
                sys.stderr.write("\r" + "P(z): computed {0:d} objects".format(i+1))
                sys.stderr.flush()

    if verbose: sys.stderr.write("\r" + "P(z): computed {0:d} objects\n".format(i+1))

    """
    write output file
    """
    if verbose: sys.stderr.write("Writing output file...")
    cols = []

    if args.key_id is not None:
        [ID], _ = getCols(fileInName[0], [args.key_id], select=fileInSelect[0])
        cols.append(fits.Column(name="ID", format='K', array=ID))

    for name in est.keys():
        cols.append(fits.Column(name=name, format='E', array=est[name]))

    #cols.append(fits.Column(name='zmedian', format='E', array=zmedian))
    cols.append(fits.Column(name='PDF',      format=str(len(PDF_ref_bins))+'E', array=pofz))

    if args.PDF_histo:
        cols.append(fits.Column(name='zTrain',      format=str(Nnei)+'E', array=zTrain))
        if len(keys_histo) > 1:
            cols.append(fits.Column(name='wTrain',      format=str(Nnei)+'E', array=wTrain))
            # cols.append(fits.Column(name='sTrain',      format=str(Nnei)+'I', array=sTrain))

    cols_bins = []
    # cols_bins.append(fits.Column(name='PDF', format=str(len(PDF_ref_bins))+'E', array=[PDF_ref_bins]))
    cols_bins.append(fits.Column(name='BINS', format='E', array=PDF_ref_bins))
    # cols_bins.append(fits.Column(name='Z_MIN', format='E', array=[PDF_ref_bins[0]]))
    # cols_bins.append(fits.Column(name='Z_MAX', format='E', array=[PDF_ref_bins[-1]]))
    # cols_bins.append(fits.Column(name='DELTA_Z', format='E', array=[PDF_ref_bins[1]-PDF_ref_bins[0]]))

    if args.merge_with_input:
        fileIn = fits.open(fileInName[0])
        hdu_0  = fileIn[0]
        if sampleSelect is not None:
            for c in fileIn[1].columns:
                cols.append(fits.Column(name=c.name, format=c.format, array=fileIn[1].data[c.name][sampleSelect]))
            hdu_1  = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        else:
            hdu_1  = fits.BinTableHDU.from_columns(fileIn[1].columns + fits.ColDefs(cols))
        if len(fileIn) > 2:
            hdu_2 = fits.BinTableHDU.from_columns(fileIn[2].columns + fits.ColDefs(cols_bins))
        else:
            hdu_2 = fits.BinTableHDU.from_columns(fits.ColDefs(cols_bins))
    else:
        hdu_0 = fits.PrimaryHDU()
        hdu_1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        hdu_2 = fits.BinTableHDU.from_columns(fits.ColDefs(cols_bins))

    hdu_1.header["EXTNAME"] = "DATA"
    hdu_2.header["EXTNAME"] = "BINS"


    if args.PDF_histo:
        hdu_1.header["z_min"] = bins[0]
        hdu_1.header["z_max"] = bins[-1]
        hdu_1.header["delta_z"] = bins[1]-bins[0]
    else:
        hdu_1.header["z_min"] = PDF_ref_bins[0]
        hdu_1.header["z_max"] = PDF_ref_bins[-1]
        hdu_1.header["delta_z"] = PDF_ref_bins[1]-PDF_ref_bins[0]


    tbhdu  = fits.HDUList([hdu_0, hdu_1, hdu_2])
    tbhdu.writeto(args.output, overwrite=True)

    if args.merge_with_input:
        fileIn.close()

    if verbose: sys.stderr.write("done\n")

    if args.plot:

        [zs], _ = getCols(fileInName[0], ["redshift"], select=fileInSelect[0])
        print "Stats (scatter, eta, bias, eta_2sig, N) =", stats(est['zmedian'], zs, [0.0, 6.0])

        # PDF_sample, PDF_sample_bins = getPDF(fileInName[0], normalise=args.no_norm, PDF_key="PDF_lephare")
        # plotPDF(pofz, PDF_ref_bins, Nsample,  "PDF.pdf",  zs=zs,  zp=est['zmedian'], PDF2=[PDF_sample, PDF_sample_bins, "lephare"])
        plotPDF(pofz, PDF_ref_bins, Nsample,  "PDF.pdf",  zs=zs,  zp=est['zmedian'])

    return





# ----------------------------------------------------- #
# Utils
# ----------------------------------------------------- #


def buildTree(RA, DEC):

    N = len(RA)

    # Spherical to euclidean geometry -> 3D
    coord3D      = np.zeros((N, 3))
    coord3D[:,0] = np.sin(np.pi/2.0 - DEC * np.pi / 180.0) * np.cos(RA * np.pi / 180.0)
    coord3D[:,1] = np.cos(np.pi/2.0 - DEC * np.pi / 180.0) * np.sin(RA * np.pi / 180.0)
    coord3D[:,2] = np.cos(np.pi/2.0 - DEC * np.pi / 180.0)


    tree = spatial.KDTree(coord3D)

    # import pickle

    # 400" in radians ~ 3 Mpc at z~1
    # Rmax = 400.0 / 3600.0 * np.pi / 180.0

    # sys.stderr.write("Searching neighbors...")
    # res = tree.query_ball_tree(tree, Rmax)
    # sys.stderr.write("done\n")

    # sys.stderr.write("Writing binary file...")
    # pickle.dump(res, open(args.output, 'w+'))
    # sys.stderr.write("done\n")

    return tree


def getCols(fileInName, keys, extension=1, selection=None, array=False, dictionary=False, interpret=True):
    '''
    Return column arrays from fileInName (fits file)

    INPUT
    fileInName: fits file name
    keys: column names list or expression (e.g. "col1 - col2")
    selection: selection (string)
    array: if set returns 2D array (N,K)
    interpret: interpret operation (default: True)

    OUTPUT
    list of column arrays
    (optional) selection array

    WARNING
    Numeric column names cannot be used
    '''

    # first open the fits file and read the data
    fileIn = fits.open(fileInName)
    data = fileIn[extension].data
    cols = fileIn[extension].columns
    fileIn.close()

    if keys is None:
        keys = cols.names

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

        if interpret:
            cmd = ""
            for s in re.split('(\W+)', k):
                if s in cols.names:
                    cmd += "data"+str_select+"['"+s+"']"
                else:
                    cmd += s

            # print cmd
            result.append(eval(cmd))
        else:
            result.append(data[k])


    if array:
        K = len(result)
        N = len(result[0])
        result_tmp = np.zeros((N,K))
        for d in range(K):
            result_tmp[:,d] = result[d]

        result = result_tmp

    if dictionary:
        result = dict(zip(keys, result))

    # return result
    if selection is not None:
        return result, select
    else:
        return result, None



def plotPDF(PDF, bins, N, fileOutName, zs=None, zp=None, PDF2=None, PDF3=None, xtitle="$z$", label="P(z)"):

    import plot_utils
    from matplotlib.backends.backend_pdf import PdfPages

    sys.stderr.write('Plot PDFs...')
    pp = PdfPages(fileOutName)

    if N > 30:
        indices = np.random.randint(0, high=N-1, size=30)
    else:
        indices = range(N)

    for i in indices:

        fig = plot_utils.plt.figure(figsize=(FIGX, FIGY)); ax = plot_utils.plt.gca(); size = FONT_SIZE

        plot_utils.axes(ax, xtitle , "PDF", [bins[0], bins[-1]], [0.0, 1.0], title=args.title)

        ylim    = max(PDF[i,:])
        nonZero = bins[PDF[i] > EPS]
        if len(nonZero) > 1:
            (xmin, xmax) = (nonZero[0], nonZero[-1])
        else:
            (xmin, xmax) = (0.0, 6.0)



        ax.fill_between(bins, 0.0, PDF[i,:], color=COLOR[0], alpha=0.5, label=label)

        if PDF2 is not None:
            ax.fill_between(PDF2[1], 0.0, PDF2[0][i,:], color=COLOR[1], label=PDF2[2], alpha=0.5)
            ylim = max(ylim, max(PDF2[0][i,:]))
            nonZero = PDF2[1][PDF2[0][i,:] > EPS]
            (xmin, xmax) = (min(nonZero[0],xmin), max(nonZero[-1],xmax))

        if PDF3 is not None:
            ax.fill_between(PDF3[1], 0.0, PDF3[0][i,:], color=COLOR[2], label=PDF3[2], alpha=0.5)
            ylim = max(ylim, max(PDF3[0][i,:]))
            nonZero = PDF3[1][PDF3[0][i,:] > EPS]
            (xmin, xmax) = (min(nonZero[0],xmin), max(nonZero[-1],xmax))

        if zp is not None:
            ax.plot([zp[i], zp[i]], [0.0, ylim], '--', color=COLOR[3], lw=2, label="$z_\mathrm{phot}$")

        if zs is not None:
            ax.plot([zs[i], zs[i]], [0.0, ylim], '--', color=COLOR[4], lw=2, label="$z_\mathrm{spec}$")

        ax.set_ylim([0.0,ylim])
        ax.set_xlim([xmin,xmax])

        ax.legend(frameon=False, loc="upper right")
        fig.set_tight_layout(True)

        pp.savefig()

    pp.close()
    sys.stderr.write('done\n')

    return


def getPDF(fileInName, normalise=False, PDF_key="PDF", select=None):
    '''
    Returns PDF array and PDF bins from input file
    '''

    thdu     = fits.open(fileInName)
    if select is not None:
        PDF      = thdu[1].data[PDF_key][select]
    else:
        PDF      = thdu[1].data[PDF_key]

    bins = thdu[2].data[PDF_key][0]
    thdu.close()
    if normalise:
        for i in range(len(PDF)):
            norm = int_trapz(bins, PDF[i,:], bins[0], bins[-1])
            if norm > EPS:
                PDF[i,:] /= norm

    return PDF, bins


def stats(zp, zs, bins):

  sigma   = np.zeros((len(bins)-1))
  eta     = np.zeros((len(bins)-1))
  eta2sig = np.zeros((len(bins)-1))
  bias    = np.zeros((len(bins)-1))
  n       = np.zeros((len(bins)-1))

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
              n[i]              = len(dist)

  return sigma, eta, bias, eta2sig, n

def max_pos_PDF(PDF_bins, PDF):
    '''
    Returns max position of PDF
    '''

    #return  PDF_bins[PDF.argmax()]

    i_peak = PDF.argmax()
    i_min  = np.maximum(0, i_peak-1)
    i_max  = i_min + 3

    if i_max > len(PDF_bins)-1:
        i_max  = np.minimum(i_peak+1, len(PDF_bins)-1)
        i_min  = i_max - 3

    new_bins  = np.arange(PDF_bins[i_min], PDF_bins[i_max], 0.001)
    inter     = interpolate.KroghInterpolator(PDF_bins[i_min:i_max], PDF[i_min:i_max])

    PDF_inter = inter(new_bins)

    return  new_bins[PDF_inter.argmax()]

def sampleFromDist(x, f, N=1, q=None):
    '''
    get random number from the distribution
    described by the function f
    f is normalised in the range between xmin and xmax
    '''

    # Whether f is a function or an array
    if hasattr(f, '__call__'):
        y = f(x)
    else :
        y = f

    # Get cumulative sum by trapeze integration
    cum_y     = np.zeros(len(x))
    cum_y[1:] = np.cumsum(np.diff(x) * (y[:-1]+y[1:])/2.0)

    # Interpolate the inverse cumulative function.
    # Division by max guarantees the cumulative function
    # varies between 0.0 and 1.0 (= normalisation)
    inv_cum_y = interpolate.interp1d(cum_y/max(cum_y), x, kind='linear')

    if q is not None:
        return inv_cum_y(q)
    else:
        if N == 1:
            return inv_cum_y(np.random.rand())
        else:
            return inv_cum_y(np.random.rand(N))

def int_trapz(x, y, a, b):
    # return the integrated value between a and b
    # use the trapeze rule
    int_range = np.linspace(a, b, num=1000)
    f         = np.interp(int_range, x, y)
    return np.sum(np.diff(int_range) * (f[:-1]+f[1:])/2.0)



def medianfromDist(x, y):

    return sampleFromDist(x, y, q=0.5)


    # Get cumulative sum by trapeze integration
    # cum_y     = np.zeros(len(x))
    # cum_y[1:] = np.cumsum(np.diff(x) * (y[:-1]+y[1:])/2.0)

    # Interpolate the inverse cumulative function.
    # Division by max guarantees the cumulative function
    # varies between 0.0 and 1.0 (= normalisation)
    # inv_cum_y = interpolate.interp1d(cum_y/max(cum_y), x, kind='linear')

    # return inv_cum_y(0.5)

def weightedSample(sample, weight):

    from scipy import stats

    # [mode, _] = stats.mode(weight)
    # weight /= mode

    indices=[]
    for i,w in enumerate(weight):
        for n in range(np.random.poisson(w)):
            indices.append(i)

    return sample[indices]


# ----------------------------------------------------- #
# Main
# ----------------------------------------------------- #

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('option',                                   help="Which quantity to plot")
    parser.add_argument('-i', '--input',  default=None,             help='input file')
    parser.add_argument('-o', '--output', default=None,             help='output file')
    parser.add_argument('-t', '--title',  type=str,   default="",   help='Title')

    parser.add_argument('-s',    '--select', type=str,   default=None, help='Selection')


    parser.add_argument('-key_id',           type=str,   default=None, help='Key of ID')

    parser.add_argument('-keys_ref',             type=str,   default=None, help='Keys to be used')
    parser.add_argument('-keys',             type=str,   default=None, help='Keys to be used')
    parser.add_argument('-keys_err',         type=str,   default=None, help='Key errors associated to keys')
    parser.add_argument('-keys_histo',       type=str,   default=None, help='Keys for building the prior [z,weight] or [z]' )

    parser.add_argument('-PDF_histo', action='store_true',          help='Set if PDF built from zspec histogram [set zspec and weight keys with -keys_histo')
    parser.add_argument('-density', action='store_true',            help='Display point density')
    parser.add_argument('-zmax',   type=float, default=4.0,         help='Upper redshift limit')

    parser.add_argument('-plot', action='store_true',               help='Plot tests')
    parser.add_argument('-output_PDF', default=None,               help='Output PDF file')

    parser.add_argument('-no_norm', action='store_false',            help='Do not normalise input (faster if input PDFs are already normalised)')

    parser.add_argument('-merge_with_input', action='store_true',            help='merge with input file if set')

    parser.add_argument('-weight_key', type=str,   default=None,     help='Weight key')


    args = parser.parse_args()

    main(args)
