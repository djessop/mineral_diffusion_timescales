#!/usr/env/python3

from scipy.stats import gamma as gamma_distn
from scipy.special import gamma, gammaln, polygamma
from scipy.optimize import fsolve, minimize, fmin

import numpy as np


class MLEGamma:
    ''' 
    Provide ML estimators for the gamma distribution, with associated
    error and other functionality.
    '''
    theta = np.array([1., 1.])
    
    def __init__(self, data):
        self.data = data

    def jacobian(self, theta, data):
        alpha, beta = theta
        x = np.array(data)
        n = len(x) 
        xbar = data.mean() 
        return np.array([n * polygamma(0, alpha) - n * np.log(beta)
                         - np.log(x).sum(),
                         n * xbar - n * alpha / beta])

    def hessian(self, theta, data):
        alpha, beta = theta
        x = np.array(data)
        n = len(x) 
        return np.array([[  n * polygamma(1, alpha), -n / beta],
                         [ -n / beta,                 n * alpha / beta**2]])

    def neg_loglikelihood(self, theta, data):
        alpha, beta = theta
        x = np.array(data)
        n = len(x)
        xbar = x.mean()
        return -(n * alpha * np.log(beta) + (alpha - 1) * np.log(x).sum() -
                 n * beta * xbar - n * gammaln(alpha))

    def fit(self, theta0, data):
        x = np.array(data)
        res = minimize(self.neg_loglikelihood,
                       x0=theta0,
                       args=(x,),
                       jac=self.jacobian)
        return res


def se_estimates(theta, x):
    '''
    Return estimates of the standard errors on parameters (shape and rate)
    of a gamma distribution as fitted to a population of data, x.

    Parameters
    ----------
    theta : list or array-like
        Fitted parameters
    x : list or array-like
        Data points

    Returns
    -------
    se_alpha, se_beta
        standard errors on alpha and beta
    '''
    alpha, beta = theta 
    n = len(x) 
    denom = np.sqrt(alpha * polygamma(1, alpha) - 1) 
    se_alpha = np.sqrt(alpha / n) / denom 
    se_beta  = np.sqrt(beta**2 * polygamma(1, alpha) / n) / denom 
    return se_alpha, se_beta 


def ml_estimates(theta0, x):
    '''
    recast as system of equations to be solved simultaneously
    f(alpha, beta) = 0
    '''
    def jac(theta, x):
        alpha, beta = theta
        n = len(x) 
        xbar = x.mean() 
        return np.array([n * np.log(beta) + np.log(x).sum() -
                         n * polygamma(0, alpha),
                         n * alpha / beta - n * xbar])

    def hes(theta, x):
        alpha, beta = theta
        n = len(x) 
        return np.array([[ -n * polygamma(1, alpha), n / beta],
                         [  n / beta, - n * alpha / beta**2]])

    
    theta = fsolve(jac, x0=(theta0), args=(x,), fprime=hes)

    return theta


