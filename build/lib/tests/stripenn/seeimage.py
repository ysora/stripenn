import argparse
import cooler
import os
import pandas as pd
import numpy as np
import warnings
import time
import sys
import cv2 as cv
import matplotlib.pyplot as plt

'''
def argumentParser():
    parser = argparse.ArgumentParser(description='Stripenn')
    parser.add_argument("cool", help="Balanced cool file.")
    parser.add_argument("position", help="Genomic position (e.g., chr1:135010000-136000000)")
    parser.add_argument("maxpixel", help="Quantile for the pixel saturation. (e.g., 0.95)")
    parser.add_argument('--norm', help="Normalization method. It should be one of the column name of Cooler.bin(). Check it with Cooler.bins().columns (e.g., KR, VC, VC_SQRT)", default='KR')
    parser.add_argument('--out', help="Output path", default='./out.png')
    args = parser.parse_args()

    cool = args.cool
    position = args.position
    out = args.out
    if out[-1] != '/':
        out += '/'
    norm = args.norm
    maxpixel = args.maxpixel
    maxpixel = float(maxpixel)

    return(cool, position, maxpixel, norm, out)
'''

def seeimage(cool, position, maxpixel, norm, out):
    #cool, position, maxpixel, norm, out = argumentParser()
    Lib = cooler.Cooler(cool)

    PossibleNorm = Lib.bins().columns
    if norm == 'None':
        norm = False
    elif norm not in PossibleNorm:
        print('Possible normalization methods are:')
        print('None')
        for n in range(3,len(PossibleNorm)):
            print(PossibleNorm[n])
        print("Invalid normalization method. Normalization method is forced to None")
        norm = False

    unbalLib = Lib.matrix(balance=norm)
    all_chromnames = Lib.chromnames
    chr = position.split(':')
    chr = chr[0]
    if chr not in all_chromnames:
        sys.exit("Invalid chromosome name.")
    cfm = unbalLib.fetch(chr)
    M = np.quantile(a=cfm[cfm > 0], q=maxpixel, interpolation='linear')
    del cfm

    A = unbalLib.fetch(position,position)

    framesize = A.shape[0]
    red = np.ones((framesize, framesize)) * 255
    blue = 255 * (M - A) / M
    blue[np.where(blue < 0)] = 0
    green = blue

    img = cv.merge((red / 255, green / 255, blue / 255))
    img = np.clip(img, a_min=0, a_max=1)


    fig = plt.subplot(111)
    plt.imshow(img)
    plt.title(position)
    plt.show()
    fig.figure.savefig(out)

    return 0
