import numpy as np
import scipy

'''
def quantile(a, prob):
    """
    Estimates the prob'th quantile of the values in a data array.

    Uses the algorithm of matlab's quantile(), namely:
    - Remove any nan values
    - Take the sorted data as the (.5/n), (1.5/n), ... , (1-.5/n) quantiles.
    - Use linear interpolation for values between (.5/n) and (1-.5/n).
    - Ues the minimum or maximum for quantiles outside that range.

    See also: scipy.stats.mstats.mquantiles
    Referred from https://stackoverflow.com/questions/13733034/equivalent-python-command-for-quantile-in-matlab
    :param a: numpy array
    :param prob: probability
    :return: quantile value
    """

    a =  np.asanyarray(a)
    a = a[np.logical_not(np.isnan(a))].ravel() # Return a contiguous flattened array
    n = a.size

    if prob >= 1 - .5/n:
        return a.max()
    elif prob <= .5 / n:
        return a.min()

    # find the two bounds we're interpreting between:
    # that is, find i such that (i+.5)/n <= prob <= (i+1.5)/n
    t = n * prob - .5
    i = int(np.floor(t))

    # partial sort so that the ith element is at position i, with bigger ones
    # to the right and smaller to the left
    a = bn.partition(a, i)

    if i == t: # did we luck out and get an interger index?
        return a[i]
    else:
        # we'll linearly interpolate between this and the next index
        smaller = a[i]
        larger = a[i+1:].min()
        if np.isinf(smaller):
            return smaller # avoid int - inf
        return smaller + (larger - smaller) * (t - i)
'''

def make_expected_vals(cooler, Ndiags):
    vals = {}
    for chrom in list(cooler.chroms()[:]['name'])[:-1]:#remove the last chrY
        mat = cooler.matrix(balance=True, sparse=True).fetch(chrom)
        diags = scipy.sparse.triu(scipy.sparse.tril(mat, Ndiags), 0, 'dia')
        for i in range(Ndiags):
            if i not in vals:
                vals[i] = []
            arr = diags.data[i][diags.data[i]==diags.data[i]]
            vals[i].append([np.mean(arr), np.std(arr), arr.shape[0]])#Std in case we want to do z-score
    avals = []
    for d in vals:
        l = np.sum([i[2] for i in vals[d]])
        m = np.array([i[0]*i[2] for i in vals[d]]).sum()/l
        s = np.array([i[1]*i[2] for i in vals[d]]).sum()/l
        avals.append((m ,s))
    return np.array(avals)

def observed_over_expected(
    matrix, mask=np.empty(shape=(0), dtype=np.bool_), dist_bin_edge_ratio=1.03
):
    """
    Normalize the contact matrix for distance-dependent contact decay.
    The diagonals of the matrix, corresponding to contacts between loci pairs
    with a fixed distance, are grouped into exponentially growing bins of
    distances; the diagonals from each bin are normalized by their average value.
    Parameters
    ----------
    matrix : np.ndarray
        A 2D symmetric matrix of contact frequencies.
    mask : np.ndarray
        A 1D or 2D mask of valid data.
        If 1D, it is interpreted as a mask of "good" bins.
        If 2D, it is interpreted as a mask of "good" pixels.
    dist_bin_edge_ratio : float
        The ratio of the largest and the shortest distance in each distance bin.
    Returns
    -------
    OE : np.ndarray
        The diagonal-normalized matrix of contact frequencies.
    dist_bins : np.ndarray
        The edges of the distance bins used to calculate average
        distance-dependent contact frequency.
    sum_pixels : np.ndarray
        The sum of contact frequencies in each distance bin.
    n_pixels : np.ndarray
        The total number of valid pixels in each distance bin.
    """
    N = matrix.shape[0]

    mask2d = np.empty(shape=(0, 0), dtype=np.bool_)
    if mask.ndim == 1:
        if mask.size > 0:
            mask2d = mask.reshape((1, -1)) * mask.reshape((-1, 1))
    elif mask.ndim == 2:
        # Numba expects mask to be a 1d array, so we need to hint
        # that it is now a 2d array
        mask2d = mask.reshape((int(np.sqrt(mask.size)), int(np.sqrt(mask.size))))
    else:
        raise ValueError("The mask must be either 1D or 2D.")

    data = np.copy(matrix).astype(np.float64)

    has_mask = mask2d.size > 0
    dist_bins = _logbins_numba(1, N, dist_bin_edge_ratio)
    dist_bins = np.concatenate((np.array([0]), dist_bins))
    n_pixels_arr = np.zeros_like(dist_bins[1:])
    sum_pixels_arr = np.zeros_like(dist_bins[1:], dtype=np.float64)
    dist_bins_max = np.max(np.where(dist_bins < 500))
    bin_idx, n_pixels, sum_pixels = 0, 0, 0

    for bin_idx, lo, hi in zip(
        range(dist_bins_max), dist_bins[:-1], dist_bins[1:]
    ):

        sum_pixels = 0
        n_pixels = 0
        for offset in range(lo, hi):
            for j in range(0, N - offset):
                if not has_mask or mask2d[offset + j, j]:
                    sum_pixels += data[offset + j, j]
                    n_pixels += 1

        n_pixels_arr[bin_idx] = n_pixels
        sum_pixels_arr[bin_idx] = sum_pixels

        if n_pixels == 0:
            continue
        mean_pixel = sum_pixels / n_pixels
        if mean_pixel == 0:
            continue

        for offset in range(lo, hi):
            for j in range(0, N - offset):
                if not has_mask or mask2d[offset + j, j]:

                    data[offset + j, j] /= mean_pixel
                    if offset > 0:
                        data[j, offset + j] /= mean_pixel

    return data, dist_bins, sum_pixels_arr, n_pixels_arr


def _logbins_numba(lo, hi, ratio=0, N=0, prepend_zero=False):
    """Make bins with edges evenly spaced in log-space.
    Parameters
    ----------
    lo, hi : int
        The span of the bins.
    ratio : float
        The target ratio between the upper and the lower edge of each bin.
        Either ratio or N must be specified.
    N : int
        The target number of bins. The resulting number of bins is not guaranteed.
        Either ratio or N must be specified.
    """
    lo = int(lo)
    hi = int(hi)
    if ratio != 0:
        if N != 0:
            raise ValueError("Please specify N or ratio")
        N = np.log(hi / lo) / np.log(ratio)
    elif N == 0:
        raise ValueError("Please specify N or ratio")
    data10 = 10 ** np.linspace(np.log10(lo), np.log10(hi), int(N))
    data10 = np.rint(data10)
    data10_int = np.sort(np.unique(data10)).astype(np.int_)
    assert data10_int[0] == lo
    assert data10_int[-1] == hi
    if prepend_zero:
        data10_int = np.concatenate((np.array([0]), data10_int))
    return data10_int

def elementwise_product_sum(K, *COL):
    ncol_K = K.shape[1]
    len_COL = len(COL)
    if ncol_K != len_COL:
        raise ValueError('The number of columns must be identical with the column size of K')
    result = []
    colSize = COL[0].size

    for i in range(1, colSize-1):
        tempVal = 0
        for j in range(0, ncol_K):
            tempVal += K[0][j]*COL[j][i-1] + K[1][j]*COL[j][i] + K[2][j]*COL[j][i+1]
        if ncol_K == 3 and tempVal < 0:
            tempVal *= -1
        result.append(tempVal)
    return result

def moving_average(x,w):
    return np.convolve(x, np.ones(w),'valid')/w
