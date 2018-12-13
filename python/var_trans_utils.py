# import all the necessary libraries
import os
import sys
import re
import pandas as pd
import numpy as np

import collections

# these are the functions to characterize and 
# plot the filters

def mean_trans(lbd, trans):
    """ Return the mean
    wavelength of the transmission
    """
    
    norm = np.trapz(trans, lbd)
    result = np.trapz(trans*lbd, lbd)/norm
    
    return result

def std_trans(lbd, trans):
    """ Return the width
    in wavelength of 
    the filter
    """

    norm = np.trapz(trans, lbd)
    mean_lambda = mean_trans(lbd, trans)
    result = np.sqrt(
        np.trapz((
            mean_lambda-lbd)**2*trans, lbd)/norm)
    
    return result

def skewness_trans(lbd, trans):
    """ Return the skewness
    in wavelength of 
    the filter
    """
        
    norm = np.trapz(trans, lbd)
    mean_lambda = mean_trans(lbd, trans)
    sigma_lambda = std_trans(lbd, trans)

    result = np.trapz(
        ((mean_lambda-lbd)/sigma_lambda)**3*trans, lbd)/norm
    
    return result

def kurtosis_trans(lbd, trans):
    """ Return the kurtosis
    in wavelength of 
    the filter
    """

    norm = np.trapz(trans, lbd)
    mean_lambda = mean_trans(lbd, trans)
    sigma_lambda = std_trans(lbd, trans)

    result = np.trapz(
        ((mean_lambda-lbd)/sigma_lambda)**4*trans, lbd)/norm
    
    return result


def moments_trans(lbd, trans):
    """ Return the first four moments
    of the transmission
    """

    result = {}
    result['mean'] = mean_trans(lbd, trans)
    result['sigma'] = std_trans(lbd, trans)
    result['skewness'] = skewness_trans(lbd, trans)
    result['kurtosis'] = kurtosis_trans(lbd, trans)
    
    return result


def plot_trans(
        df, title, xlim, ylim=(0.0,1.0), ax=None):
    """ plot filters in 
    a data frame and return 
    the figure object 
    """
    
    ax = df.plot.line(
        x='lambda', stacked=False, legend=False,
        title=title, xlim=xlim,
        ylim=ylim, ax=ax)

    ax.set_xlabel(r'$\lambda$ [\AA{}]')
    ax.set_ylabel(r'Transmission')

    return ax

def plot_stats(
        df_center, dfoff_center, roff_center, lbd, 
        title, xlabel='r [mm]', 
        stats=['mean', 'sigma', 'skewness', 'kurtosis'], 
        ax=None, ylim=None, legend_size = 'small', 
        plot_legend = False, file_out = None):
    """ plot the filter stats compared to 
    the center
    """
    
    # profile in the center averaged over the 7 positions
    trans_center = df_center.mean(axis=1)
    
    # stats at the center
    mean_lambda_center = mean_trans(lbd, trans_center)
    sigma_lambda_center = std_trans(lbd, trans_center)
    skewness_lambda_center = skewness_trans(lbd, trans_center)
    kurtosis_lambda_center = kurtosis_trans(lbd, trans_center)
    
    mean_lambda = []
    sigma_lambda = []
    skewness_lambda = []
    kurtosis_lambda = []

    for o in dfoff_center.keys():
        trans = dfoff_center[o]
        mean_lambda.append(mean_trans(lbd, trans))
        sigma_lambda.append(std_trans(lbd, trans))
        skewness_lambda.append(skewness_trans(lbd, trans))
        kurtosis_lambda.append(kurtosis_trans(lbd, trans))
    
    if ax is None:
        fig, ax = plt.subplots()
        
    if 'mean' in stats:
        ax.scatter(
            roff_center, mean_lambda-mean_lambda_center, s=100, 
           label='Mean', marker='o')
    if 'sigma' in stats:
        ax.scatter(
            roff_center, sigma_lambda-sigma_lambda_center, s=100, 
            label='Sigma', marker='s')
    if 'skewness' in stats:
        ax.scatter(
            roff_center, 100.0*(
                skewness_lambda-skewness_lambda_center), s=100, 
            label=r'Skewness$\times100$', marker = '^')

    if 'kurtosis' in stats:
        ax.scatter(
            roff_center, 100.0*(
                kurtosis_lambda-kurtosis_lambda_center), s=100, 
            label=r'Kurtosis$\times100$', marker = '^')
        
    if file_out is not None:
        stats_dict = collections.OrderedDict({})
        stats_dict['mean'] = mean_lambda-mean_lambda_center
        stats_dict['sigma'] = sigma_lambda-sigma_lambda_center
        stats_dict['skewness'] = skewness_lambda-skewness_lambda_center
        stats_dict['kurtosis'] = kurtosis_lambda-kurtosis_lambda_center
        df = pd.DataFrame(stats_dict)
        df.to_csv(file_out, index=False)
        

    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\Delta\lambda$ [\AA{}]')
    
    if plot_legend:
        ax.legend(loc='upper left', prop={'size': legend_size}, frameon=True)   
        

    ax.plot(ax.get_xlim(), [0.0, 0.0], color='black', linestyle=':')

    # ax.legend(bbox_to_anchor=(1.1, 1.1))
    
    return ax




    # these are the functions to create and 
# change the transmission shape

def create_r_top_hat_trans(name, top=0.5):
    """ Create a top-hat transmission 
    simlar to a r-band """
    N = 10000
    
    result = {}
    if name is 'r' :
        result['lambda'] = np.linspace(5000.0, 8000.0, N)
        result['trans'] = np.zeros(N)
        result['trans'][(5600.0 < result['lambda']) \
                        & (result['lambda'] < 7000.0)] = top
    return result

def shift_trans(lbd, trans, shift):
    """ Shift the filter
    transmission in wavelength.
    A positive shift goes towards
    higher wavelengths (red) 
    
    INPUT
    - lbd: wavelength
    - trans: transmission values
    - shift: in wavelength

    OUTPUT
    - modified filter
    """
    
    return np.interp(lbd-shift, lbd, trans)

def stretch_trans(lbd, trans, stretch):
    """ Modify the transmission by 
    stretching the transmission 
    around the mean.
    
    The strecth value is the fractional 
    increase of the transmission,
    between -1 and 1.
    
    INPUT
    - lbd: wavelength
    - trans: transmission values
    - stretch (between -1 and 1)

    OUTPUT
    - modified filter
    """
    
    result = np.zeros(len(trans))
    
    # record where the transmission is non zero
    pos = np.argwhere(trans > 0.1).flatten()
    
    # mean lambda
    mean = mean_trans(lbd, trans)
    
    lbd_stretched = (lbd-mean)*(1.0+stretch)+mean
    
    result = np.interp(lbd, lbd_stretched, trans)
    
    return result

def tilt_trans(lbd, trans, tilt):
    """ Modify the transmission by 
    tilting the top around the mean.
    
    The tilt value is the fractional 
    increase of the transmission at
    the red edge, between -0.5 and 0.5.
    
    A positive tilt will skew 
    the filter to the red.
    
    INPUT
    - lbd: wavelength
    - trans: transmission values
    - tilt (between -0.5 and 0.5)

    OUTPUT
    - modified filter
    """
    
    result = np.zeros(len(trans))
    
    # record where the transmission is non zero
    pos = np.argwhere(trans > 0.1).flatten()
    
    # mean lambda
    mean = mean_trans(lbd, trans)
    
    # sigma lambda
    sigma = std_trans(lbd, trans)

    # fractional increase
    frac_increase = 1.0+tilt/sigma \
        * (lbd[pos]-mean)

    # compute transmission
    result[pos] = trans[pos]*frac_increase

    return result

def softening_trans(lbd, trans, softening):
    """ Modify the transmission by 
    softening the edges through 
    a convolution with a trapezoid
    
    The softening value is the width of the 
    trapezoid relative to the sigma of 
    the distribution
    
    INPUT
    - lbd: wavelength
    - trans: transmission values
    - softening (between 0 and 1)

    OUTPUT
    - modified filter
    """
    
    # result = np.zeros(len(trans))

    if softening == 0.0:
        return trans
    
    # mean lambda
    mean = mean_trans(lbd, trans)

    # sigma lambda
    sigma = std_trans(lbd, trans)
    
    # kernel = trapezoid
    # x = np.array([-1.0, -0.1, +0.1,  +1.0])*sigma*softening
    # y = np.array([0.0, 1.0, 1.0, 0.0])
    x = np.array([-1.0, 0.0,  +1.0])*sigma*softening
    y = np.array([0.0, 1.0, 0.0])


    kernel_x = np.arange(x[0], x[-1], lbd[1]-lbd[0])
    kernel_y = np.interp(kernel_x, x, y)
    
    result = np.convolve(trans, kernel_y, 'same')
    
    return result/np.max(result)*np.max(trans)