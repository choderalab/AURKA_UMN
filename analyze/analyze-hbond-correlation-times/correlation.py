#=============================================================================================
# IMPORTS
#=============================================================================================
import math
import numpy as np
import numpy.linalg
from pymbar.utils import ParameterError

#=============================================================================================
# Issue warning on import.
#=============================================================================================

LongWarning = "Warning on use of the timeseries module: If the inherent timescales of the system are long compared to those being analyzed, this statistical inefficiency may be an underestimate.  The estimate presumes the use of many statistically independent samples.  Tests should be performed to assess whether this condition is satisfied.   Be cautious in the interpretation of the data."

#sys.stderr.write(LongWarning + '\n')

#=============================================================================================
# METHODS
#=============================================================================================

#=============================================================================================

def integrate_autocorrelation_function(C_n, tvec):
    """Integrate a normalized fluctuation autocorrelation function to get an integrated correlation time.
    """

    mintime = 3

    # Accumulate the integrated correlation time by computing the normalized correlation time at
    # increasing values of t.  Stop accumulating if the correlation function goes negative, since
    # this is unlikely to occur unless the correlation function has decayed to the point where it
    # is dominated by noise and indistinguishable from zero.

    T = tvec.max()
    N = len(C_n)
    g = 1.0
    for n in range(N):
        C = C_n[n]
        t = tvec[n]

        # Terminate if the correlation function has crossed zero and we've computed the correlation
        # function at least out to 'mintime'.
        if (C <= 0.0) and (n > mintime):
            break

        if n == 0:
            increment = tvec[0]
        else:
            increment = tvec[n] - tvec[n-1]

        # Accumulate contribution to the statistical inefficiency.        
        g += 2.0 * C * (1.0 - float(t) / float(T)) * float(increment)

    # g must be at least unity
    if (g < 1.0):
        g = 1.0

    # g = 1 + 2*tau
    # tau = (g - 1)/2
    tau = (g - 1.0) / 2.0

    # Return the computed correlation time tau.
    return tau
#=============================================================================================


def statisticalInefficiencyMultiple(A_kn, fast=False, return_correlation_function=False):
    """Estimate the statistical inefficiency from multiple stationary timeseries (of potentially differing lengths).

    Parameters
    ----------
    A_kn : list of np.ndarrays
        A_kn[k] is the kth timeseries, and A_kn[k][n] is nth value of timeseries k.  Length is deduced from arrays.

    fast : bool, optional, default=False
        f True, will use faster (but less accurate) method to estimate correlation
        time, described in Ref. [1] (default: False)
    return_correlation_function : bool, optional, default=False
        if True, will also return estimates of normalized fluctuation correlation function that were computed (default: False)

    Returns
    -------
    g : np.ndarray,
        g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
        We enforce g >= 1.0.
    Ct : list (of tuples)
        Ct[n] = (t, C) with time t and normalized correlation function estimate C is returned as well if return_correlation_function is set to True

    Notes
    -----
    The autocorrelation of the timeseries is used to compute the statistical inefficiency.
    The normalized fluctuation autocorrelation function is computed by averaging the unnormalized raw correlation functions.
    The fast method described in Ref [1] is used to compute g.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
        histogram analysis method for the analysis of simulated and parallel tempering simulations.
        JCTC 3(1):26-41, 2007.

    Examples
    --------

    Estimate statistical efficiency from multiple timeseries of different lengths.

    >>> from pymbar import testsystems
    >>> N_k = [1000, 2000, 3000, 4000, 5000]
    >>> tau = 5.0 # exponential relaxation time
    >>> A_kn = [ testsystems.correlated_timeseries_example(N=N, tau=tau) for N in N_k ]
    >>> g = statisticalInefficiencyMultiple(A_kn)

    Also return the values of the normalized fluctuation autocorrelation function that were computed.

    >>> [g, Ct] = statisticalInefficiencyMultiple(A_kn, return_correlation_function=True)

    """

    # Convert A_kn into a list of arrays if it is not in this form already.
    if (type(A_kn) == np.ndarray):
        A_kn_list = list()
        if A_kn.ndim == 1:
            A_kn_list.append(A_kn.copy())
        else:
            [K, N] = A_kn.shape
            for k in range(K):
                A_kn_list.append(A_kn[k, :].copy())
        A_kn = A_kn_list

    # Determine number of timeseries.
    K = len(A_kn)

    # Get the length of each timeseries.
    N_k = np.zeros([K], np.int32)
    for k in range(K):
        N_k[k] = A_kn[k].size

    # Compute average timeseries length.
    Navg = np.array(N_k, np.float64).mean()

    # Determine total number of samples.
    N = np.sum(N_k)

    # Initialize statistical inefficiency estimate with uncorrelated value.
    g = 1.0

    # Compute sample mean.
    mu = 0.0
    for k in range(K):
        mu += np.sum(A_kn[k])
    mu /= float(N)

    # Construct and store fluctuation timeseries.
    dA_kn = list()
    for k in range(K):
        dA_n = A_kn[k] - mu
        dA_kn.append(dA_n.copy())

    # Compute sample variance from mean of squared fluctuations, to ensure that C(0) = 1.
    sigma2 = 0.0
    for k in range(K):
        sigma2 += np.sum(dA_kn[k] ** 2)
    sigma2 /= float(N)

    # Initialize statistical inefficiency estimate with uncorrelated value.
    g = 1.0

    # Initialize storage for correlation function.
    Ct = list()  # Ct[n] is a tuple (t, C) of the time lag t and estimate of normalized fluctuation correlation function C

    # Accumulate the integrated correlation time by computing the normalized correlation time at
    # increasing values of t.  Stop accumulating if the correlation function goes negative, since
    # this is unlikely to occur unless the correlation function has decayed to the point where it
    # is dominated by noise and indistinguishable from zero.
    t = 1
    increment = 1
    while (t < N_k.max() - 1):
        # compute unnormalized correlation function
        numerator = 0.0
        denominator = 0.0
        for k in range(K):
            if (t >= N_k[k]):
                continue  # skip trajectory if lag time t is greater than its length
            dA_n = dA_kn[k]  # retrieve trajectory
            x = dA_n[0:(N_k[k] - t)] * dA_n[t:N_k[k]]
            numerator += np.sum(x)  # accumulate contribution from trajectory k
            denominator += float(x.size)  # count how many overlapping time segments we've included

        C = numerator / denominator

        # compute normalized fluctuation correlation function at time t
        C = C / sigma2
        # print "C[%5d] = %16f (%16f / %16f)" % (t, C, numerator, denominator)

        # Store estimate of correlation function.
        Ct.append((t, C))

        # Terminate if the correlation function has crossed zero.
        # Note that we've added a hack (t > 10) condition to avoid terminating too early in correlation functions that have a strong negative peak at
        if (C <= 0.0) and (t > 10):
            break

        # Accumulate contribution to the statistical inefficiency.
        g += 2.0 * C * (1.0 - float(t) / Navg) * float(increment)

        # Increment t and the amount by which we increment t.
        t += increment

        # Increase the interval if "fast mode" is on.
        if fast:
            increment += 1

    # g must be at least unity
    if (g < 1.0):
        g = 1.0

    # Return statistical inefficency and correlation function estimate, if requested.
    if return_correlation_function:
        return (g, Ct)

    # Return the computed statistical inefficiency.
    return g
    
#=============================================================================================

def unnormalizedFluctuationCorrelationFunction(A_n, B_n=None, N_max=None, dot_product_function=None):
    """Compute the unnormalized fluctuation (cross) correlation function of (two) stationary timeseries using a specified dot product function

    C(t) = <A(0).B(t)>

    This may be useful in diagnosing odd time-correlations in timeseries data.

    Parameters
    ----------
    A_n : np.ndarray
        A_n[n] is nth value of timeseries A.  Length is deduced from vector.
    B_n : np.ndarray
        B_n[n] is nth value of timeseries B.  Length is deduced from vector.
    N_max : int, default=None
        if specified, will only compute correlation function out to time lag of N_max
    dot_product_function : function, optional, default=None
        if specified, this function will be used to evaluate A_t0 . B_t1

    Returns
    -------
    C_n : np.ndarray
        C_n[n] is the normalized fluctuation auto- or cross-correlation function for timeseries A(t) and B(t).

    Notes
    -----
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    This procedure may be slow.
    The statistical error in C_n[n] will grow with increasing n.  No effort is made here to estimate the uncertainty.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

    """

    # If B_n is not specified, set it to be identical to A_n.
    if B_n is None:
        B_n = A_n

    # Get the length of the timeseries.
    N = A_n.size

    # Set maximum time to compute correlation functon for.
    if (not N_max) or (N_max > N - 1):
        N_max = N - 1

    # Be sure A_n and B_n have the same dimensions.
    if(A_n.shape != B_n.shape):
        raise ParameterError('A_n and B_n must have same dimensions.')

    # allocate storage for normalized fluctuation correlation function
    C_n = np.zeros([N_max + 1], np.float64)

    # Compute unnormalized correlation function.
    if dot_product_function:
        for t in range(0, N_max + 1):
            C_n[t] = 0.0
            denominator = 0.0
            for t0 in range(0, N - t):
                C_n[t] += dot_product_function(A_n[t0], B_n[t0+t])
                C_n[t] += dot_product_function(B_n[t0], A_n[t0+t])
                denominator += 2.0
            C_n[t] /= denominator
    else:
        raise Exception("dot_product_function must be specified")

    # Return the computed correlation function
    return C_n
#=============================================================================================


def unnormalizedFluctuationCorrelationFunctionMultiple(A_kn, B_kn=None, N_max=None, dot_product_function=None, nskip=1):
    """Compute the unnormalized fluctuation (cross) correlation function of (two) timeseries from multiple timeseries samples.

    C(t) = <A(0).B()t>
    This may be useful in diagnosing odd time-correlations in timeseries data.

    Parameters
    ----------
    A_kn : Python list of numpy arrays
        A_kn[k] is the kth timeseries, and A_kn[k][n] is nth value of timeseries k.  Length is deduced from arrays.
    B_kn : Python list of numpy arrays
        B_kn[k] is the kth timeseries, and B_kn[k][n] is nth value of timeseries k.  B_kn[k] must have same length as A_kn[k]
    N_max : int, optional, default=None
        if specified, will only compute correlation function out to time lag of N_max
    dot_product_function : function, optional, default=None
        if specified, this function will be used to evaluate A_t0 . B_t1
    nskip : int, optional, default=1
        Compute C(t) at this interval in frames

    Returns
    -------
    tvec : np.ndarray of int
        Frame intervals for which C_t is evaluated
    C_n[n] : np.ndarray
        The unnormalized auto- or cross-correlation function for timeseries A(t) and B(t).

    Notes
    -----
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    This procedure may be slow.
    The statistical error in C_n[n] will grow with increasing n.  No effort is made here to estimate the uncertainty.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

    """

    # If B_kn is not specified, define it to be identical with A_kn.
    if B_kn is None:
        B_kn = A_kn

    # Ensure the same number of timeseries are stored in A_kn and B_kn.
    if (len(A_kn) != len(B_kn)):
        raise ParameterError("A_kn and B_kn must contain corresponding timeseries -- different numbers of timeseries detected in each.")

    # Determine number of timeseries stored.
    K = len(A_kn)

    # Ensure both observable trajectories in each timeseries are of the same length.
    for k in range(K):
        A_n = A_kn[k]
        B_n = B_kn[k]
        if A_n.size != B_n.size:
            raise ParameterError("A_kn and B_kn must contain corresponding timeseries -- lack of correspondence in timeseries lenghts detected.")

    # Get the length of each timeseries.
    N_k = np.zeros([K], np.int32)
    for k in range(K):
        N_k[k] = A_kn[k].size

    # Determine total number of samples.
    N = np.sum(N_k)

    # Set maximum time to compute correlation functon for.
    if (not N_max) or (N_max > max(N_k) - 1):
        N_max = max(N_k) - 1

    # Compute lag times for which C(t) will be evaluated
    tvec = np.arange(0, N_max+1, nskip)
    nval = len(tvec)

    # allocate storage for normalized fluctuation correlation function
    C_n = np.zeros([nval], np.float64)
    N_n = np.zeros([nval], np.float64)

    # Compute unnormalized correlation function.
    if dot_product_function:
        for n in range(nval):
            t = tvec[n] # lag time at which computation is performed
            C_n[n] = 0.0
            N_n[n] = 0.0
            for k in range(K):
                for t0 in range(0, N_k[k] - t, nskip):
                    try:
                        C_n[n] += 0.5 * dot_product_function(A_kn[k][t0], B_kn[k][t0+t])
                    except Exception as e:
                        print(k)
                        print(t0)
                        print(t)
                        print(k,t0)
                        print(k,t0+t)
                        raise(e)
                    C_n[n] += 0.5 * dot_product_function(B_kn[k,t0], A_kn[k,t0+t])
                    N_n[n] += 1.0
            C_n[n] /= N_n[n]
    else:
        raise Exception("dot_product_function must be specified")

    # Return the computed correlation function
    return [tvec, C_n, N_n]

#========================
