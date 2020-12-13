import numpy as np
import cv2 as cv
from scipy.signal import convolve2d
import math
import scipy.ndimage
import time

def imBrightness3D(img, In=([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), Out=([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])):
    # "J = low_out +(high_out - low_out).* ((I - low_in)/(high_in - low_in)).^ gamma"
    # Modified from this code: https://www.programmersought.com/article/32635116380/

    if img.max() > 1 or img.min() < 0:
        raise ValueError('Pixel values must be rescaled to zero to one.')

    def imgconvert(simg, h, w, k, low_in, high_in, low_out, high_out):
        imgOut = np.zeros((h,w))
        [x_low, y_low] = np.where(simg <= low_in)
        [x_high, y_high] = np.where(simg > high_in)
        [x_mid, y_mid] = np.where((simg > low_in) & (simg <= high_in))
        imgOut[x_low, y_low] = low_out
        imgOut[x_high, y_high] = high_out
        imgOut[x_mid, y_mid] = k * (simg[x_mid,y_mid] - low_in) + low_out
        return imgOut

    ([r_low_in, g_low_in, b_low_in], [r_high_in, g_high_in, b_high_in]) = In
    ([r_low_out, g_low_out, b_low_out], [r_high_out, g_high_out, b_high_out]) = Out


    r_k = (r_high_out - r_low_out) / (r_high_in - r_low_in)
    g_k = (g_high_out - g_low_out) / (g_high_in - g_low_in)
    b_k = (b_high_out - b_low_out) / (b_high_in - b_low_in)

    h, w = img.shape[:2]

    r_imgOut = imgconvert(img[:,:,0], h, w, r_k, r_low_in, r_high_in, r_low_out, r_high_out)
    g_imgOut = imgconvert(img[:,:,1], h, w, g_k, g_low_in, g_high_in, g_low_out, g_high_out)
    b_imgOut = imgconvert(img[:,:,2], h, w, b_k, b_low_in, b_high_in, b_low_out, b_high_out)

    imgOut = cv.merge((r_imgOut, g_imgOut, b_imgOut))

    return imgOut

def auto_canny(image, sigma = 0.33):
    """
    Canny edge detection without lower- and upper-bound setting.
    Referred from: https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
    :param image: gray scale image
    :param sigma: I don't know
    :return: edge image
    """
    # compute the median of the single channel pixel intensities
    v = np.median(image)

    # apply automatic Canny edge detection using the computed median
    lower = int(max(0, (1.0 - sigma) * v))
    upper = int(min(255, (1.0 + sigma) * v))
    edged  = cv.Canny(image, lower, upper)

    return edged

def verticalLine(M, L=60, H=120):
    """
    :param M: A binary matrix consisting of 0 and 1.
    :return: Pixel coordinates with vertical line
    """
    Gx = [[-1, 0, 1],[-2, 0, 2],[-1, 0, 1]]
    Gy = [[1, 2, 1],[0, 0, 0],[-1, -2, -1]]

    Filtered_X = convolve2d(M, Gx, mode='same')
    Filtered_Y = convolve2d(M, Gy, mode='same')

    orientation = np.arctan2(Filtered_X, Filtered_Y)
    orientation = orientation * 180/math.pi
    [x_neg, y_neg] = np.where(orientation < 0)
    orientation[x_neg, y_neg] += 360
    [x_one, y_one] = np.where((orientation > L) & (orientation < H))
    #[_, y_two] = np.where((orientation > L+180) & (orientation < H+180))
    y_one = np.subtract(y_one,1)
    del orientation
    orientation = np.zeros(M.shape)
    orientation[x_one, y_one] = 1
    #orientation[x_two, y_two] = 1
    return orientation
'''
    # Direction and orientations
    orientation = np.zeros(M.shape)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ori = math.atan2(Filtered_Y[i,j],Filtered_X[i,j])
            ori = ori * 180 / math.pi
            if (ori >= H*(-1) and ori <= L*(-1)) or (ori >= L and ori <= H):
                orientation[i,j] = 1
            else:
                orientation[i,j] = 0

    orientation2 = np.arctan2(Filtered_X, Filtered_Y)
    sum(sum(orientation == orientation))
'''

def block(mat,c):
    """

    :param mat: 2D array
    :param c: column index
    :return: t, END
    """

    count = 0;
    MAX = 0;
    S = mat.shape[1]
    #value = mat[:,c]

    #if c > 0:
    #    value_L = mat[:, c-1]
    #else:
    #    value_L = np.zeros(S)
    #if c < S-1:
    #    value_R = mat[:, c+1]
    #else:
    #    value_R = np.zeros(S)

    value = mat[:,max(0,c-1):min(S,c+2)]
    value = value.max(axis = 1)
    END = 0
    #END_REV = 0
    buffer = 0
    for i in range(S):
        if value[i] == 1:
            count += 1
            J = i
        elif buffer < 5:
            buffer += 1
        else:
            if count > MAX:
                MAX = count
                END = J
            count = 0
            buffer = 0
    if count > MAX:
        MAX = count
        END = J
    '''
    for i in range(S):
        if value[i,:].max() == 1:
            count += 1
            J = i
        elif buffer < 5:
            buffer += 1
        else:
            if count > MAX:
                MAX = count
                END = J
            count = 0
            buffer = 0

    if count > MAX:
        MAX = count
    t = MAX
    '''
    '''
    t0 = time.time()
    count_rev = 0
    MAX_REV = 0
    buffer = 0
    for i in range(S-1, -1, -1):
        if value[i,:].max() == 1:
            count_rev += 1
            J = i
        elif buffer < 5:
            buffer += 1
        else:
            if count_rev > MAX_REV:
                MAX_REV = count_rev
                END_REV = J
            count_rev = 0
            buffer = 0

    if count_rev > MAX_REV:
        MAX_REV = count_rev
        END_REV = J
    print(time.time()-t0)
    test1 = abs(c - END)
    test2 = abs(c - END_REV)
    
    if test1 > test2:
        t = MAX
        END = END
    else:
        t = MAX_REV
        END = END_REV
    '''
    t = MAX
    if END < c:
        END = END - t + 1
    return t, END

def Canny(img, threshold):

    imgsize = img.shape
    nrow = imgsize[0]
    ncol = imgsize[1]

    # Magic numbers
    PercentOfPixelsNotEdges = .7 # Used for selecting thresholds
    ThresholdRatio = .4          # Low thresh is this freaction of the high.
    thresh = threshold

    # Calculate gradients using a derivative of Gaussian
    dx, dy = smoothGradient(img)

    # Calculate Magnitude of Gradient
    magGrad = hypot(dx, dy)

    # Normalize for threshold selection
    magmax = magGrad.max()
    if magmax > 0:
        magGrad = magGrad / magmax

    # Calculate directions/orientation
    arah = np.zeros(imgsize)
    vert = np.zeros(imgsize)
    arah2 = np.zeros(imgsize)
    for i in range(nrow):
        for j in range(ncol):
            ori = math.atan2(dy[i,j],dx[i,j])
            ori = ori * 180 / math.pi
            if ori < 0:
                ori += 360
            arah[i,j] = ori

    for i in range(nrow):
        for j in range(ncol):
            if (arah[i,j] >= 0 and arah[i,j] < 22.5) or (arah[i,j] >=157.5 and arah[i,j] < 202.5) or (arah[i,j] >= 337.5 and arah[i,j] <=360):
                arah2[i,j] = 0
            elif (arah[i,j] >= 22.5 and arah[i,j] < 67.5) or (arah[i,j] >= 202.5 and arah[i,j] < 247.5):
                arah2[i,j] = 45
            elif (arah[i,j] >=67.5 and arah[i,j] < 112.5) or (arah[i,j] >= 247.5 and arah[i,j] < 292.5):
                arah2[i,j] = 90
            elif (arah[i,j] >= 112.5 and arah[i,j] < 157.5) or (arah[i,j] >= 292.5 and arah[i,j] < 337.5):
                arah2[i,j] = 135

    BW = np.zeros(imgsize)
    for i in range(1,nrow-1,1):
        for j in range(1,ncol-1,1):
            if arah2[i,j] == 0:
                BW[i,j] = magGrad[i,j] == max(magGrad[i,j], magGrad[i,j+1], magGrad[i,j-1])
            elif arah2[i,j] == 45:
                BW[i,j] = magGrad[i,j] == max(magGrad[i,j], magGrad[i+1,j-1], magGrad[i-1,j+1])
            elif arah2[i,j] == 90:
                BW[i,j] = magGrad[i,j] == max(magGrad[i,j], magGrad[i+1,j], magGrad[i-1,j])
            elif arah2[i,j] == 135:
                BW[i,j] = magGrad[i,j] == max(magGrad[i,j], magGrad[i+1, j+1], magGrad[i-1,j-1])

    BW = np.multiply(BW, magGrad)

    # Hysteresis Thresholding
    T_Low = thresh * ThresholdRatio * BW.max()
    T_High = thresh * BW.max()
    T_res = np.zeros(imgsize)

    for i in range(nrow):
        for j in range(ncol):
            if BW[i,j] < T_Low:
                T_res[i,j] = 0
            elif BW[i,j] > T_High:
                T_res[i,j] = 1
            elif BW[i+1,j]>T_High or BW[i-1,j]>T_High or BW[i,j+1]>T_High or BW[i,j-1]>T_High or BW[i-1, j-1]>T_High or BW[i-1, j+1]>T_High or BW[i+1, j+1]>T_High or BW[i+1, j-1]>T_High:
                T_res[i,j] = 1

    return T_res



def hypot(a, b):
    if type(a) != np.ndarray:
        a = np.array(a)
    if type(b) != np.ndarray:
        b = np.array(b)
    if a.shape != b.shape:
        raise ValueError('Two arrays should have same dimension!')
    res = a**2 + b**2
    res = res ** 0.5
    return res

def smoothGradient(I, sigma=math.sqrt(2)):
    """

    :param I: Image object
    :param sigma: Standard deviation of the filter, specified as a numeric scalar. Default = sqrt(2)
    :return: dx, dy
    """
    # Determine filter length
    filterExtent = math.ceil(4 * sigma)
    x = list(range(-1*filterExtent, filterExtent+1, 1))

    # Create 1-D Gaussian Kernel
    c = 1/(math.sqrt(2*math.pi)*sigma)
    gaussKernel = [c * math.exp(-(i**2)/(2*sigma**2)) for i in x]

    # Normalize to ensure kernel sums to one
    gaussKernel = [i/sum(gaussKernel) for i in gaussKernel]

    # Create 1-D Derivative of Gauss Kernel
    derivGaussKernel = simple_gradient(gaussKernel)

    # Normalize to ensure kernel sums to zero
    negVals = derivGaussKernel < 0
    posVals = derivGaussKernel > 0
    derivGaussKernel[posVals] = derivGaussKernel[posVals]/sum(derivGaussKernel[posVals])
    derivGaussKernel[negVals] = derivGaussKernel[negVals]/abs(sum(derivGaussKernel[negVals]))

    gaussKernel = np.array([gaussKernel])
    derivGaussKernel = np.array([derivGaussKernel])
    # Compute smoothed numerical gradient of image I along x (horizontal)
    # deriction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
    # version of image I.
    GX = scipy.ndimage.convolve(I, np.transpose(gaussKernel), mode='nearest')
    GX = scipy.ndimage.convolve(GX, derivGaussKernel, mode='nearest')

    # Compute smoothed numerical gradient of image I along y (vertical)
    # direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
    # version of image I.
    GY = scipy.ndimage.convolve(I, gaussKernel, mode='nearest')
    GY = scipy.ndimage.convolve(GY, np.transpose(derivGaussKernel), mode='nearest')

    return GX, GY

def simple_gradient(f):
    rowflag = False
    ndim = 1
    indx = [len(f), 1]
    loc = [[i for i in range(indx[0])],[1]]
    siz = indx

    # first dimension
    g = np.zeros(siz[0])
    h = loc[0]
    n = siz[0]

    # take forward differences on left and right edges
    if n > 1:
        g[0] = (f[1] - f[0]) / (h[1] - h[0])
        g[n-1] = (f[n-1] - f[n-2])/(h[n-1] - h[n-2])

    if n > 2:
        g[1:n-1] = [(f[i]-f[j])/(h[k]-h[l]) for i,j,k,l in zip(range(2,n), range(0,n-2), range(2,n), range(0,n-2))]

    return g

