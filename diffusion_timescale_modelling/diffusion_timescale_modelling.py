#!/usr/bin/env python
# coding: utf-8


from functools import partial

from scipy.special import erfc, gammaln, polygamma, gamma
from scipy.stats import expon, kstest, norm, chi2 
from scipy.stats import gamma as gamma_distn
from scipy.optimize import curve_fit, newton, fmin
from scipy.constants import R  # Ideal-gas constant
from matplotlib.ticker import MultipleLocator

from ml_stats.mle_gamma import se_estimates


import pandas as pd # Use pandas' inherrent excel support, rather than numpy
import numpy as np
import matplotlib.pyplot as plt

from xlrd import XLRDError
import xlrd


pd.options.display.float_format = '{:,.2f}'.format


fO2 = 6.31e-7    # Partial pressure of O2 for NNO buffer with dNNO=0.8, see
                 # Pichavant et al., (2018).
E   = -308000    # Activation energy
T   = 1281.75    # Melt matrix temperature, from Dohmen et al. (2016)
#D0  =  2.8551e-7 # Base diffusion coefficient


def opx_model(fO2=fO2):
    '''
    Return D0, the diffusion coefficient (at absolute zero), for the opx 
    mineral system using the model from Dohmen et al, (2016).

    Parameters
    ----------
    fO2 : float
        partial pressure of O2 in Pa.  Default value is 6.31e-7 Pa, valid 
        for dNNO = 0.8, T ~ 900°C as per Pichavant et al. (2018).

    Returns
    -------
    D0 : float
        Diffusion coefficient (at absolute zero), as required by the Arhenius 
        equation.  If x is in Pa, D0 is returned in m2/s.

    SEE ALSO
    --------
    diffusion_coeff
    temp_from_diff_coeff
    '''
    return 1.12e-6 * fO2 ** 0.053


def diffusion_coeff(T, fO2=fO2):
    '''
    Return the temperature-dependent diffusion coefficient according to the 
    current model.  Note that the activation energy, E, is defined as being 
    negative.

    Parameters
    ----------
    T : float
        Absolute temperature/[K]
    fO2 : float
        Partial pressure of O2/[Pa]

    SEE ALSO
    --------
    opx_model
    temp_from_diff_coeff
    '''
    D0 = opx_model(fO2)
    
    return D0 * np.exp(E / (R * T))


def temp_from_diff_coeff(D, fO2=fO2):
    '''
    Return temperature from diffusion coefficient according to the 
    current model, the inverse of the diffusion_coeff function.  Note that
    E is defined as being -ve.

    SEE ALSO
    --------
    diffusion_coeff
    opx_model
    '''
    D0 = opx_model(fO2)

    return E / (R * np.log(D/D0))


def analytical_soln(x, C1, C2, loc, tau, D): 
    '''
    Model for elemental diffusion, following Girona and Costa (2013)
    
    Accounts for the upper and lower concentration limits (C1, C2), the 
    midpoint of the curve (loc), the timescale (tau) and (optionally) 
    the diffusion coefficient.  
    
    Specify which mode of fit (adjust for diffusivity or not) through the 
    length of the list of initial guesses that is passed to __curve_fit__.
    '''
    
    return C2 + (C1 - C2) / 2 * erfc((x - loc) / (2 * np.sqrt(D * tau)))


def fit_wrapper(x, y, p0, bounds, jac=None, function_to_fit=analytical_soln):
    # Length of upper and lower bounds have to be identical
    lower_bounds, upper_bounds = bounds
    assert len(lower_bounds) == len(upper_bounds)

    if len(p0) == 5:
        # The analytical solution is non linear so we need a good guess at 
        # the solution in order for the solver to converge.
        popt, pcov = curve_fit(function_to_fit, x, y, p0, 
                               bounds=(lower_bounds, upper_bounds))
        diffusion  = popt[4]
        Test       = temp_from_diff_coeff(popt[4])
    else:
        # Ensure that list of initial values has no entry for diffusion - this 
        # will cause the solver to ignore diffusion (and hence temperature) as
        # ajustable parameters.
        assert len(p0) == 4
        if len(upper_bounds) == 5:
            upper_bounds.pop()
            lower_bounds.pop()
        D = diffusion_coeff(T, fO2)
        fn_wrapper = partial(function_to_fit, D=D)
        popt, pcov = curve_fit(fn_wrapper, x, y, p0, jac=jac,
                               bounds=(lower_bounds, upper_bounds))
        diffusion  = D
        Test       = T
        
    timescale = popt[3] / (3600 * 24)  # In days
    ts_sigma  = 1.96 * np.sqrt(np.diag(pcov))[3] / (3600 * 24)
    result_tuple = (timescale, ts_sigma, Test, diffusion)

    return popt, pcov, result_tuple


def coefficient_determination(xdata, ydata, popt, model=analytical_soln):
    '''
    Calculate the "Pearson" regression coefficient of determination, R2
    '''
    residuals = ydata - model(xdata, *popt)
    ybar      = ydata.mean()
    SSres     = (residuals**2).sum()
    SStot     = ((ydata - ybar)**2).sum()

    return 1 - SSres/SStot


def eval_ts_data(df, R2thresh=.85, chi2rthresh=False, n_sigma=1., 
                path_to_params='./model_fits',
                path_to_data='./SEM_traverse_data/', data_fmt='xls'):
    '''
    Examine data from SEM traverses and disqualify any that do not satisfy 
    goodness-of-fit and model accuracy tests.

    Parameters
    ----------
    df : pandas.DataFrame
    R2thresh : scalar
    chi2thresh : boolean
    n_sigma : scalar
    path_to_params : string
        location of optimised parameter files
    path_to_data : string
        location of data files, using glob format
    data_fmt : string
        format of data files, e.g. xls, xlsx, csv, dat.

    Returns
    -------
    keep_list : array-like
        list of data to keep
    
    '''
    from glob import glob

    def abs_normalise(x, data):
        '''
        return the absolute normalised-standard value of a variable given 
        a distribution of data
        '''
        return np.abs((x - data.mean()) / data.std())

    keep_list = []
    for row in df.itertuples():
        filename = ('./SEM_traverse_data/' + row.eruption + '/' +
                    row.sample + '.xls')
        data     = pd.read_excel(filename, sheet_name='raw')
        x, y     = data['distance'], data['greyscale']
        popt     = np.load('./model_fits/' + row.sample + '-popt.npy')

        loc      = popt[2]
        # Eliminate data if 'loc' far from centre (i.e. 'shoulders' not 
        # captured), or if the R2 value is too low.  Komo's curvature test?
        if (abs_normalise(loc, x) > 1.):
            keep_list.append(False)
        elif row.R2 < R2thresh:
            keep_list.append(False)
        else:
            keep_list.append(True)

    return np.array(keep_list)


def unique_sample_names(df):

    samples = df['sample']
    cropped_samples = []
    for sample in samples:
        cropped_samples.append(sample[:-2])

    return sorted(list(set(cropped_samples)))


def sorted_data_to_df(df, R2thresh=0.85):
    keep_list = eval_ts_data(df, R2thresh=R2thresh)
    df      = df[keep_list]
    
    samples = unique_sample_names(df)

    sort_df = pd.DataFrame()

    for sample in samples:
        sub_df = df[df['sample'].str.contains(sample)]
        d = sub_df.mean().to_dict()
        d.pop('R2')
        d.pop('ts_sigma')
        d['n_good']   = len(sub_df)  # only count "good" R2 values
        d['ts_std']   = sub_df.std()['timescale']
        d['eruption'] = sub_df['eruption'].iloc[0]
        d['sample']   = sample
        #d['R2'] = sub_df['R2'].mean()
        sort_df = pd.concat([sort_df,
                             pd.DataFrame.from_dict(d, orient='index').T])
      
    return sort_df.sort_values(by=['eruption', 'timescale']).reset_index(
        drop=True)[['eruption', 'sample', 'Tmeas', 'Test',
                    'timescale', 'n_good', 'ts_std']]


def run_model_fitting(filenames, do_plot=True,
                      sheetname = 'Dan, WH37 processing, usabl (2)'):
    from xlrd import XLRDError
    import xlrd

    # A list for storing dicts.  This will later be converted to a dataframe
    # and used for analysing and averaging the timescales.
    list_of_samples = []

    fname_error_log = open('filename_error.log', 'w')
    fit_error_log   = open('fit_error.log', 'w')
    key_error_log   = open('key_error.log', 'w')
    ran_files_log   = open('ran_files.log', 'w')
    
    for filename in filenames:
        try:
            data = pd.read_excel(filename, sheet_name='raw')
            x    = data['distance']
            y    = data['greyscale']
            ws   = xlrd.open_workbook(filename).sheet_by_name(sheetname)
            T    = ws.cell(3, 13).value

            D    = diffusion_coeff(T)
            Dlow = diffusion_coeff(T - 30)
            Dupp = diffusion_coeff(T + 30)
    
            # Bounds for optimisation are: 
            lower_bounds = [-np.inf, -np.inf, -np.inf, 0, Dlow]
            upper_bounds = [ np.inf,  np.inf,  np.inf, np.inf, Dupp]
            p0   = [y.max(), y.min(), x.mean(), 2e6, D]

    
            code, eruption, *foo = filename.split('/')[::-1]
            code = code[:-4]

            ##--------------------DATA FITTING--------------------##
            popt, pcov, fitted_parameters = fit_wrapper(
                x, y, p0, bounds=[lower_bounds, upper_bounds])
            timescale, ts_sigma, Test, diff_coef = fitted_parameters

            np.save('./model_fits/%s-popt.npy' % code, popt)
            np.save('./model_fits/%s-pcov.npy' % code, pcov)

            if do_plot:
                data, R2, figax = plot_data_model(filename)
            else:
                data = pd.read_excel(filename, sheet_name='raw')
                R2   = coefficient_determination(data[:,0], data[:,1], popt)
            # Update the list of dicts with information for each sample.
            list_of_samples.append({'eruption': eruption, 'sample': code,
                                    'Tmeas': T, 'Test': Test,
                                    'timescale': timescale,
                                    'ts_sigma': ts_sigma,
                                    'R2': R2})

        except(XLRDError):
            mesg = 'File error with %s\n' % filename
            print(mesg, end='')
            fname_error_log.write(mesg)
        except(KeyError):
            mesg = 'Key not found in file %s\n' % filename
            print(mesg, end='')
            key_error_log.write(mesg)
        except(RuntimeError):
            mesg = 'Could not fit case %s\n' % code
            print(mesg, end='')
            fit_error_log.write(mesg)

        ran_files_log.write(filename + '\n')

    df = pd.DataFrame(list_of_samples)
    df.to_csv('timescale_fitting_df.csv', index=False, float_format='%.4f')

    fname_error_log.close()
    fit_error_log.close()   
    key_error_log.close()
    ran_files_log.close()

    return df


def plot_data_model(filename, popt=None, pcov=None, savefig=True,
                    sheetname='Dan, WH37 processing, usabl (2)'):
    '''
    
    '''

    code, eruption, *foo = filename.split('/')[::-1]
    code = code[:-4]

    data = pd.read_excel(filename, sheet_name='raw')
    x, y = data['distance'], data['greyscale']

    all_data = read_raw_data(filename)

    if popt is None:
        popt = np.load('./model_fits/%s-popt.npy' % code)
    if pcov is None:
        pcov = np.load('./model_fits/%s-pcov.npy' % code)

    if len(popt) == 4:
        print(filename)
        print(popt)

    timescale = popt[3] / (3600 * 24)  # in days
    ts_sigma  = np.sqrt(np.diag(pcov)[3]) / (3600 * 24)

    diff_coef = popt[4]
    Test      = temp_from_diff_coeff(diff_coef)


    # Create a figure for plotting
    fig, ax = plt.subplots()
    ax.plot(all_data[:,0] * 1e6, all_data[:,1], '.', label='all data')
    ax.plot(x * 1e6, y, '.r', label='selected data')
                
    _ = ax.set_xlabel(r'Traverse location/[$\mu$m]', fontsize=14)
    _ = ax.set_ylabel(r'Greyscale intensity', fontsize=14)

    model_soln = analytical_soln(x, *popt).values
    ax.plot(x * 1e6, model_soln, '-k', label='model fit')
    ax.legend(loc='right')

    # The midpoint of the curve
    ax.axvline(popt[2] * 1e6, c='gray', zorder=1)

    # Fill the span of mean +/- std of data
    ax.axvspan((x.mean() - x.std()) * 1e6,
               (x.mean() + x.std()) * 1e6,
               color='grey', zorder=0, alpha=.5)
    ax.axvline((x.mean() - x.std()) * 1e6, ls='--', c='grey',
               zorder=0, alpha=.6)
    ax.axvline((x.mean() + x.std()) * 1e6, ls='--', c='grey',
               zorder=0, alpha=.6)

    timescale_text = ('Timescale = %.2f $\pm$ %.2f days\nD = ' \
                      '%g m$^2$/s\nTemperature = %.2f K' % 
                      (timescale, ts_sigma, diff_coef, Test))
    TSkwargs = dict(x=.05, y=.05,
                    horizontalalignment='left',
                    verticalalignment='bottom')
    GFkwargs = dict(x=.95, y=.95,
                    horizontalalignment='right',
                    verticalalignment='top')

    # Goodness of fit.  N.B. Reduced Chi-squared should be close to 
    # unity for a good fit but ONLY for weighted Chi-squared statistic.
    misfit  = ((analytical_soln(x, *popt) - y)**2).sum()  
    DF      = len(data) - len(popt)
    chi2est = misfit / DF
    R2      = coefficient_determination(x, y, popt)

    alpha = 0.6
    if model_soln[0] < model_soln[-1]:
        TSkwargs.update(y=.95, verticalalignment='top')
        GFkwargs.update(y=.05, verticalalignment='bottom')
        
    ax.text(s=timescale_text, **TSkwargs, 
            transform = ax.transAxes,
            bbox=dict(boxstyle='round', alpha=alpha, color='w'))
        
    ax.text(s='Misfit = %.2f\n$R^2$ = %.4f\n' % (misfit, R2), 
            transform = ax.transAxes, **GFkwargs,
            bbox=dict(boxstyle='round', alpha=alpha, color='w'))
        
    ax.set_title(code)
        
    if savefig:
        fig.savefig(filename[:-4] + '.png')
        plt.close()
        
    return data, R2, (fig, ax)


def calc_fit_plot(filename='./22-06-sorted-timescales.xlsx',
                  sheet_name='1657CE', timescalekey='timescale',
                  ppf_head=.01, ppf_tail=.999, ppf_cutoff=None,
                  ts_cutoff=None, cdf=False, nbins=10, figsize=(8, 6)):
    '''
    Fit and plot a gamma distribution to data
    '''
    # Read and trim the data, if required
    data = pd.read_excel(filename, 
                         sheet_name=sheet_name)[timescalekey].values.ravel()
    if ts_cutoff is not None:
        data = data[data <= ts_cutoff]

    # Alpha value for transparency of annotation boxes
    alpha = .5
    
    # Open a figure and axes
    fig, ax = plt.subplots(figsize=figsize)

    # Create a histogram over which we plot the predicted distribution...
    out = ax.hist(data, density=True, cumulative=cdf, bins=nbins, 
                  histtype='bar', fc='#1f77b4a0', ec='#1f77b4ff')
    # ...then fit the data to the distribution using MLE 
    popt = gamma_distn.fit(data)

    # Get limits of the histogram...
    xlims = list(ax.get_xlim())
    ylims = list(ax.get_ylim())
    # ...and calculate some offsets as a percentage of the axes' spans
    xoffs = np.diff(xlims) * .01
    yoffs = np.diff(ylims) * .02
    
    # Plot the theoretical distribution over the histogram
    x = np.linspace(gamma_distn.ppf(ppf_head, *popt),
                    gamma_distn.ppf(ppf_tail, *popt), 101)
    if cdf:
        distn = gamma_distn.cdf
        yl    = 'cdf'
        scale = 1.05
    else:
        distn = gamma_distn.pdf
        yl    = 'pdf'
        scale = 1.25
    label = 'gamma ' + yl
    ax.plot(x, distn(x, *popt), 'r-', lw=2, label=label)
    
    # Shade area between distribution curve and abscissa, up to the ppf cutoff
    if ppf_cutoff is not None:
        x0, x1 = (gamma_distn.ppf(ppf_head, *popt),
                  gamma_distn.ppf(ppf_cutoff, *popt))
        xx = np.linspace(x0, x1, 101)
        yf = distn(xx, *popt)
        ax.plot([x1]*2, [0, yf[-1]], ls='--', c='r', lw=2, 
                alpha=alpha)
        ax.annotate(ppf_cutoff, xy=(x1, yf[-1]*scale), color='r',
                    fontweight='bold', ha='left', va='bottom',
                    bbox=dict(boxstyle='round', alpha=alpha, color='w'))
        ax.fill_between(xx, yf, color='r', alpha=alpha/2)

    
    # Expected value of the pdf (N.B. popt = (shape, loc, scale))
    distn_mean   = gamma_distn.mean(*popt)
    distn_median = gamma_distn.median(*popt)
    yc_mean      = distn(distn_mean, *popt)
    yc_median    = distn(distn_median, *popt)
    
    # Indicate and annotate key features of the CDF
    ax2_order = 0
    if cdf:
        ax.plot([xlims[0], distn_median], [.5]*2, c='r', ls='--', lw=2)
        ax.plot([distn_median]*2, [0, yc_median], '--r', lw=2)
        ax.annotate('$\mathrm{median}(x) = %.1f$ d' % distn_median, 
                    xy=(distn_median*scale, yc_mean/2), color='r',
                    fontweight='normal', ha='left', va='center', 
                    bbox=dict(boxstyle='round', alpha=alpha, color='w'),
                    rotation=90, zorder=5)
        xlims = ax.set_xlim(xlims)
        # Also add the pdf of the gamma distn on secondary axis
        ax2 = ax.twinx()
        ax2.plot(x, gamma_distn.pdf(x, *popt), '-g', lw=2, zorder=1)
        ax2.set_ylabel('PDF', fontsize=16, color='g')
        ax2_order = ax2.get_zorder()

    
    # Standard error of expectation
    alpha, beta = [popt[0], 1/popt[2]]
    se_alpha, se_beta = se_estimates([alpha, beta], data)
    se_tau = np.sqrt((se_alpha/alpha)**2 + (se_beta/beta)**2)

    ax.plot([distn_mean]*2, [0, yc_mean], '-r', lw=3)
    ax.plot([xlims[0], distn_mean], [yc_mean]*2, '-r', lw=3)
    ax.annotate('$E[x] = %.1f\pm %.3f$ d' % (distn_mean, se_tau),
                xy=(distn_mean*scale, yc_mean/2),
                color='r', fontweight='normal',
                ha='left', va='center',
                bbox=dict(boxstyle='round', alpha=alpha, color='w'),
                rotation=90, zorder=5+ax2_order)
    ax.annotate('$%.2f$' % yc_mean, xy=(xlims[0]+xoffs, yc_mean+yoffs), 
                fontweight='normal', ha='left', va='bottom', color='r',
                bbox=dict(boxstyle='round', alpha=alpha, color='w'),
                zorder=5+ax2_order)
    
    # Indicate the spreadsheet name on the plot
    ax.annotate(sheet_name, xy=(.18, .90), xycoords='axes fraction', 
                bbox=dict(boxstyle='round', alpha=alpha, color='w'),
                fontweight='bold', ha='right', zorder=2)
    
    # Perform the Kolmogorov-Smirnov test for goodness-of-fit
    KS_result = kstest(data, 'gamma', args=popt)
    ax.annotate('pvalue = %.3f' % KS_result[1], xy=(.9, .1), 
                xycoords='axes fraction', 
                fontweight='bold', ha='right', va='bottom',
                bbox=dict(boxstyle='round', alpha=alpha, color='w'),
                zorder=5+ax2_order)
    
    # Plot labels and prettification
    ax.set_xlabel(r'$\tau$/[days]', fontsize=16)
    ax.set_ylabel(yl.upper(), fontsize=16, color='r')
    ax.yaxis.set_minor_locator(MultipleLocator(.1))

    fig.savefig('timescale_distributions_%s.png' % sheet_name)
    
    return popt, data, (se_alpha, se_beta, se_tau), (fig, ax)


def read_raw_data(filename, sheetname='Dan, WH37 processing, usabl (2)',
                  cols=(3,4), start_row=8):

    ws = xlrd.open_workbook(filename).sheet_by_name(sheetname)
    data = np.array([ws.col_values(cols[0])[start_row:],
                     ws.col_values(cols[1])[start_row:]]).T

    return data


if __name__ == '__main__':
    from glob import glob

    filenames = [f for f in sorted(glob('SEM_traverse_data/**',
        recursive=True)) if f.endswith('.xls')]
    # filenames = [f for f in sorted(glob('SEM_traverse_data/720BCE/**',
    #     recursive=True)) if f.endswith('.xls')]

    df = run_model_fitting(filenames)
