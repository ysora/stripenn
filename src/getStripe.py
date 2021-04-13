import numpy as np
import pandas as pd
import math
import matplotlib
import cv2 as cv
from src import stats
from src import ImageProcessing
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import statistics as stat
from skimage import feature
import time
import random
from scipy import signal
from joblib import Parallel, delayed
from tqdm import tqdm

class getStripe:
    def __init__(self, unbalLib, resol, minH, maxW, canny, all_chromnames, chromnames, all_chromsizes, chromsizes, core):
        self.unbalLib = unbalLib
        self.resol = resol
        self.minH = minH
        self.maxW = maxW
        self.canny = canny
        self.all_chromnames = all_chromnames
        self.all_chromsizes = all_chromsizes
        self.chromnames = chromnames
        self.chromsizes = chromsizes
        self.core = core
        self.chromnames2sizes={}
        for i in range(len(self.all_chromnames)):
            self.chromnames2sizes[self.all_chromnames[i]] = self.all_chromsizes[i]

    def ExpectedCount(self, coolinfo, ChrList):
        res = {}
        chrom_size = coolinfo.chromsizes
        chrom_cum_size = np.nancumsum(chrom_size)
        chrom_names = chrom_size.keys()
        nbin = coolinfo.binsize
        chridx = [c for c in range(len(chrom_names)) if chrom_names[c] in ChrList]

        for ci in chridx:
            CHROM = chrom_names[ci]
            SIZE = chrom_size[ci]
            if ci == 0:
                r_start = 0
            else:
                r_start = chrom_cum_size[ci - 1] + 1
            r_start = int(np.ceil(r_start / nbin))
            r_end = r_start + int(np.ceil(SIZE / nbin))

            mean_list = []

            for i in range(0,401):
                val = np.mean(np.diag(self.unbalLib[r_start:r_end,r_start:r_end], i))
                mean_list.append(val)
            res[CHROM] = mean_list
        return res

    def getQuantile(self, coolinfo, ChrList, quantile):
        def SplitVal(w):
            with np.errstate(divide='ignore', invalid='ignore'):
                w_start = 5000 * w
                w_start += r_start
                w_end = 5000 * w + 4999
                if w_end >= np.floor(CHROMSIZE / nbin):
                    w_end = int(np.floor(CHROMSIZE / nbin))
                w_end += r_start

                w_start = int(w_start)
                w_end = int(w_end)
                temp = self.unbalLib[w_start:w_end, w_start:w_end][self.unbalLib[w_start:w_end, w_start:w_end] > 0]
                return temp

        res = {}
        chrom_size = coolinfo.chromsizes
        chrom_cum_size = np.nancumsum(chrom_size)
        chrom_names = chrom_size.keys()
        nbin = coolinfo.binsize
        # Index for selected chromosomes
        chridx = [c for c in range(len(chrom_names)) if chrom_names[c] in ChrList]
        chridx.sort()
        for ci in chridx:
            CHROM = chrom_names[ci]
            CHROMSIZE = chrom_size[ci]
            if ci == 0:
                r_start = 0
            else:
                r_start = chrom_cum_size[ci - 1]

            r_start = np.ceil(r_start / nbin) + 1

            hs = np.empty(0)

            N_windows = int(np.ceil(CHROMSIZE / nbin / 5000))

            W = Parallel(n_jobs=1)(delayed(SplitVal)(ww) for ww in tqdm(range(N_windows)))
            for i in range(len(W)):
                hs = np.concatenate((hs, W[i]))
            qt = np.quantile(hs, quantile)
            res[CHROM] = qt
            del(hs)
        return res

    def getQuantile_slow(self, coolinfo, ChrList, quantile):
        def SplitVal(k,w,r_start,r_end):
            with np.errstate(divide='ignore', invalid='ignore'):
                w_start = k * w
                w_start += r_start
                w_end = (k + 1) * w - 1
                if w_end >= np.floor(CHROMSIZE / nbin):
                    w_end = int(np.floor(CHROMSIZE / nbin))
                w_end += r_start

                w_start = int(w_start)
                w_end = int(w_end)
                r_start = int(r_start)
                r_end = int(r_end)
                temp = self.unbalLib[w_start:w_end, r_start:r_end][self.unbalLib[w_start:w_end, r_start:r_end] > 0]
                return temp

        res = {}
        chrom_size = coolinfo.chromsizes
        chrom_cum_size = np.nancumsum(chrom_size)
        chrom_names = chrom_size.keys()
        nbin = coolinfo.binsize
        # Index for selected chromosomes
        chridx = [c for c in range(len(chrom_names)) if chrom_names[c] in ChrList]
        chridx.sort()
        for ci in chridx:
            CHROM = chrom_names[ci]
            CHROMSIZE = chrom_size[ci]
            L = int(np.ceil(CHROMSIZE / nbin))
            if ci == 0:
                r_start = 0
            else:
                r_start = chrom_cum_size[ci - 1]

            r_end = chrom_cum_size[ci]

            r_start = np.ceil(r_start / nbin) + 1
            r_end = np.ceil(r_end / nbin)

            hs = np.empty(0)

            w = int(np.floor(25000000 / L))

            N_windows = int(np.ceil(L / w))

            W = Parallel(n_jobs=1)(delayed(SplitVal)(ww,w,r_start,r_end) for ww in tqdm(range(N_windows)))
            for i in range(len(W)):
                hs = np.concatenate((hs, W[i]))
            qt = np.quantile(hs, quantile)
            res[CHROM] = qt
            del(hs)
        return res

    def getQuantile_original(self, coolinfo, ChrList, quantile):

        res = {}
        chrom_size = coolinfo.chromsizes
        chrom_names = chrom_size.keys()
        # Index for selected chromosomes
        chridx = [c for c in range(len(chrom_names)) if chrom_names[c] in ChrList]

        chridx.sort()
        for ci in chridx:
            CHROM = chrom_names[ci]

            mat = self.unbalLib.fetch(CHROM)
            qt = np.quantile(mat[mat>0], quantile)
            res[CHROM] = qt
            del mat
        return res

    def mpmean(self):

        def calc(n):
            pixels_mean = [0 for x in range(framesize)]
            counts_mean = [0 for x in range(framesize)]
            start = self.resol * n * framesize + 1
            end = self.resol * (n + 1) * framesize
            end2 = self.resol * (n + 2) * framesize
            if end > chrsize:
                end = chrsize
            if end2 > chrsize:
                end2 = chrsize
            rows = str(chr) + ":" + str(start) + "-" + str(end)
            cols = str(chr) + ":" + str(start) + "-" + str(end2)
            cfm = self.unbalLib.fetch(rows, cols)
            cfm_rows = cfm.shape[0]
            cfm_cols = cfm.shape[1]
            cfm_max = np.max(cfm)

            for i in range(cfm_rows):
                for j in range(min(framesize, cfm_cols)):
                    if i + j >= cfm_cols:
                        continue
                    else:
                        count = cfm[i, i + j]
                        if np.isnan(count):
                            count = 0
                        pixels_mean[j] += count
                        counts_mean[j] += 1
            del cfm
            return pixels_mean, counts_mean

       #meantable = [[] for x in range(len(self.chromnames))]
        meantable = {}

        for chridx in range(len(self.chromnames)):
            chr = self.chromnames[chridx]
            print('Processing Chromosome: ' + str(chr))
            chrsize = self.chromnames2sizes[chr]
            rowsize = int(np.ceil(chrsize / self.resol))
            framesize = 400
            nframes = math.ceil(rowsize / framesize)

            with np.errstate(divide='ignore'):
                result = Parallel(n_jobs=self.core)(delayed(calc)(n) for n in tqdm(range(nframes)))

            means = []
            for i in range(len(result[0][1])):
                pixelsum = [result[j][0][i] for j in range(nframes)]
                countsum = [result[j][1][i] for j in range(nframes)]
                pixelsum = sum(pixelsum)
                countsum = sum(countsum)
                meanval = pixelsum/countsum
                means.append(meanval)
            meantable[chr] = means

        return meantable


    def nulldist(self):
        with np.errstate(divide='ignore', invalid='ignore'):
            t_background_start = time.time()

            samplesize = (self.all_chromsizes / np.sum(self.all_chromsizes)) * 1000
            samplesize = np.uint64(samplesize)
            dif = 1000 - np.sum(samplesize)
            notzero = np.where(samplesize != 0)
            chromnames2 = [self.all_chromnames[i] for i in notzero[0]]
            samplesize[0] += dif

            def available_cols(chr):
                with np.errstate(divide='ignore', invalid='ignore'):
                    chr = str(chr)
                    chrsize = self.chromnames2sizes[chr]
                    itera = min(chrsize/self.resol/500, 25)
                    itera = np.uint64(itera)
                    unitsize = np.floor(chrsize/self.resol/itera)
                    unitsize = np.uint64(unitsize)
                    poolsum = 0
                    for it in range(itera):
                        test_region_start = np.uint64(unitsize * self.resol * it+1)
                        test_region_end = np.uint64(unitsize * self.resol * (it+1))
                        if test_region_start > test_region_end:
                            a = test_region_start
                            test_region_start = test_region_end
                            test_region_end = a
                        position = str(chr) + ":" + str(test_region_start) + "-" + str(test_region_end)
                        mat = self.unbalLib.fetch(position, position)
                        mat = nantozero(mat)
                        #mat = np.round(mat, 1)
                        matsum = np.sum(mat, axis=1)
                        zeroindex = np.where(matsum == 0)
                        poolsum += (len(matsum) - len(zeroindex[0]))

                return poolsum

            print('1. Calculating the number of available columns ... \n')
            n_available_col = Parallel(n_jobs=self.core)(delayed(available_cols)(chr) for chr in tqdm(chromnames2))

            samplesize = (n_available_col/ np.sum(n_available_col)) * 1000
            samplesize = np.uint64(samplesize)
            dif = 1000 - np.sum(samplesize)
            notzero = np.where(samplesize != 0)
            chromnames2 = [chromnames2[i] for i in notzero[0]]
            samplesize[0] += dif

            def main_null_calc(chr):
                with np.errstate(divide='ignore', invalid='ignore'):
                    background_size = 50000/self.resol
                    background_up = np.floor(background_size / 2)
                    background_down = background_size - background_up
                    background_up = int(background_up)
                    background_down = int(background_down)
                    background_size = int(background_size)
                    chr = str(chr)
                    #c = np.where(chromnames2 == chr)[0]
                    c = chromnames2.index(chr)

                    # Modified in Dec 11 2020
                    ss = samplesize[c]
                    chrsize = self.chromnames2sizes[chr]
                    itera = min(chrsize/self.resol/500, 25)
                    itera = int(itera)
                    unitsize = np.floor(chrsize/self.resol/itera)
                    unitsize = np.uint64(unitsize)

                    bgleft_up = np.zeros((400, 0))
                    bgright_up = np.zeros((400, 0))
                    bgleft_down = np.zeros((400, 0))
                    bgright_down = np.zeros((400, 0))

                    n_pool = []

                    for it in range(itera):
                        sss = int(ss/itera)
                        test_region_start1 = int(unitsize * self.resol * it+1)
                        test_region_start0 = test_region_start1 - (400*self.resol)
                        test_region_end1 = int(unitsize * self.resol * (it+1))
                        test_region_end2 = int(unitsize * self.resol * (it+1)+(400*self.resol))
                        if test_region_end2 > chrsize:
                            test_region_end2 = chrsize-1
                        if test_region_end1 > chrsize - 400*self.resol:
                            test_region_end1 = chrsize - 400*self.resol
                        if test_region_start0 <= 1:
                            test_region_start0 = 1
                        test_region_start0 = int(test_region_start0)

                        position1 = str(chr) + ":" + str(test_region_start1) + "-" + str(test_region_end1)
                        position2 = str(chr) + ":" + str(test_region_start0) + "-" + str(test_region_end2)
                        mat = self.unbalLib.fetch(position1, position2)
                        mat = nantozero(mat)
                       # mat = np.round(mat, 1)
                        nrow = mat.shape[0]
                        matsum = np.sum(mat, axis=1)
                        zeroindex = np.where(matsum == 0)
                        pool = [x for x in list(range(nrow)) if x not in zeroindex[0].tolist()]
                        pool = [x for x in pool if x > 20 and x < (unitsize - 20)]
                        n_pool.append(len(pool))
                        if len(pool) == 0:
                            del mat
                        elif len(pool) < sss:
                            randval = random.choices(pool, k=len(pool))
                            tableft_up = np.zeros((400, len(pool)))
                            tabcenter_up = np.zeros((400, len(pool)))
                            tabright_up = np.zeros((400, len(pool)))
                            tableft_down = np.zeros((400, len(pool)))
                            tabcenter_down = np.zeros((400, len(pool)))
                            tabright_down = np.zeros((400, len(pool)))

                            for i in range(len(pool)):
                                x = randval[i]
                                for j in range(0,400):
                                    #det = np.random.choice([0,1])
                                    y_down = x + j
                                    y_up = x - j

                                    #if det == 0:
                                    #    y = x + j
                                    #else:
                                    #    y = x - j
                                    if it > 0 :
                                        y_down = y_down + 400
                                        y_up = y_up + 400

                                    tableft_up[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_up - background_up):(y_up + background_down)])
                                    tabcenter_up[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_up - background_up):(y_up + background_down)])
                                    tabright_up[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_up - background_up):(y_up + background_down)])

                                    tableft_down[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_down - background_up):(y_down + background_down)])
                                    tabcenter_down[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_down - background_up):(y_down + background_down)])
                                    tabright_down[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_down - background_up):(y_down + background_down)])

                            bgleft_up_temp = np.subtract(tabcenter_up, tableft_up)
                            bgright_up_temp = np.subtract(tabcenter_up, tabright_up)
                            bgleft_up = np.column_stack((bgleft_up, bgleft_up_temp))
                            bgright_up = np.column_stack((bgright_up, bgright_up_temp))
                            bgleft_down_temp = np.subtract(tabcenter_down, tableft_down)
                            bgright_down_temp = np.subtract(tabcenter_down, tabright_down)
                            bgleft_down = np.column_stack((bgleft_down, bgleft_down_temp))
                            bgright_down = np.column_stack((bgright_down, bgright_down_temp))

                            del mat
                        else:
                            randval = random.choices(pool, k=sss)
                            tableft_up = np.zeros((400, sss))
                            tabcenter_up = np.zeros((400, sss))
                            tabright_up = np.zeros((400, sss))
                            tableft_down = np.zeros((400, sss))
                            tabcenter_down = np.zeros((400, sss))
                            tabright_down = np.zeros((400, sss))

                            for i in range(sss):
                                x = randval[i]
                                for j in range(0, 400):
                                    y_down = x + j
                                    y_up = x - j
                                   # det = np.random.choice([0,1])
                                   # if det == 0:
                                   #     y = x + j
                                   # else:
                                   #     y = x - j
                                    if it > 0:
                                        y_up = y_up + 400
                                        y_down = y_down + 400
                                    tableft_up[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_up - background_up):(y_up + background_down)])
                                    tabcenter_up[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_up - background_up):(y_up + background_down)])
                                    tabright_up[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_up - background_up):(y_up + background_down)])
                                    tableft_down[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_down - background_up):(y_down + background_down)])
                                    tabcenter_down[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_down - background_up):(y_down + background_down)])
                                    tabright_down[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_down - background_up):(y_down + background_down)])

                            bgleft_up_temp = np.subtract(tabcenter_up, tableft_up)
                            bgright_up_temp = np.subtract(tabcenter_up, tabright_up)
                            bgleft_down_temp = np.subtract(tabcenter_down, tableft_down)
                            bgright_down_temp = np.subtract(tabcenter_down, tabright_down)

                            del mat
                            bgleft_up = np.column_stack((bgleft_up, bgleft_up_temp))
                            bgright_up = np.column_stack((bgright_up, bgright_up_temp))
                            bgleft_down = np.column_stack((bgleft_down, bgleft_down_temp))
                            bgright_down = np.column_stack((bgright_down, bgright_down_temp))

                    depl = int(ss - bgleft_up.shape[1])
                    if depl > 0:
                        rich = np.argmax(n_pool)
                        test_region_start1 = int(unitsize * self.resol * rich+1)
                        test_region_start0 = int(test_region_start1 - (400*self.resol))
                        test_region_end1 = int(unitsize * self.resol * (rich+1))
                        test_region_start2 = int(unitsize * self.resol * rich+1)
                        test_region_end2 = int(unitsize * self.resol * (rich+1)+(400*self.resol))
                        if test_region_end2 > chrsize:
                            test_region_end2 = chrsize-1
                        if test_region_end1 > chrsize - 400*self.resol:
                            test_region_end1 = chrsize - 400*self.resol
                        if test_region_start0 <= 1:
                            test_region_start0 = 1

                        position1 = str(chr) + ":" + str(test_region_start1) + "-" + str(test_region_end1)
                        position2 = str(chr) + ":" + str(test_region_start0) + "-" + str(test_region_end2)
                        mat = self.unbalLib.fetch(position1, position2)
                        mat = nantozero(mat)
                        #mat = np.round(mat, 1)
                        nrow = mat.shape[0]
                        matsum = np.sum(mat, axis=1)
                        zeroindex = np.where(matsum == 0)
                        pool = [x for x in list(range(nrow)) if x not in zeroindex[0].tolist()]
                        pool = [x for x in pool if x > 20 and x < (unitsize - 20)]
                        randval = random.choices(pool, k=depl)
                        tableft_up = np.zeros((400, depl))
                        tabcenter_up = np.zeros((400, depl))
                        tabright_up = np.zeros((400, depl))
                        tableft_down = np.zeros((400, depl))
                        tabcenter_down = np.zeros((400, depl))
                        tabright_down = np.zeros((400, depl))
                        for i in range(depl):
                            x = randval[i]
                            for j in range(0, 400):
                                y_down = x + j
                                y_up= x - j
                                #det = np.random.choice([0, 1])

                                #if det == 0:
                                #    y = x + j
                                #else:
                                #    y = x - j
                                if it > 0:
                                    y_up = y_up + 400
                                    y_down = y_down+400
                                    tableft_up[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_up - background_up):(y_up + background_down)])
                                    tabcenter_up[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_up - background_up):(y_up + background_down)])
                                    tabright_up[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_up - background_up):(y_up + background_down)])
                                    tableft_down[j, i] = np.mean(mat[(x - background_up - background_size):(x - background_up), (y_down - background_up):(y_down + background_down)])
                                    tabcenter_down[j, i] = np.mean(mat[(x - background_up):(x + background_down), (y_down - background_up):(y_down + background_down)])
                                    tabright_down[j, i] = np.mean(mat[(x + background_down):(x + background_down + background_size), (y_down - background_up):(y_down + background_down)])

                        bgleft_up_temp = np.subtract(tabcenter_up, tableft_up)
                        bgright_up_temp = np.subtract(tabcenter_up, tabright_up)
                        bgleft_down_temp = np.subtract(tabcenter_down, tableft_down)
                        bgright_down_temp = np.subtract(tabcenter_down, tabright_down)

                        del mat
                        bgleft_up = np.column_stack((bgleft_up, bgleft_up_temp))
                        bgright_up = np.column_stack((bgright_up, bgright_up_temp))
                        bgleft_down = np.column_stack((bgleft_down, bgleft_down_temp))
                        bgright_down = np.column_stack((bgright_down, bgright_down_temp))

                        return bgleft_up, bgright_up, bgleft_down, bgright_down
            # apply parallel.
            print('2. Constituting background ... \n')
            result = Parallel(n_jobs=self.core)(delayed(main_null_calc)(chr) for chr in tqdm(chromnames2))
            bgleft_up = np.zeros((400,0))
            bgright_up = np.zeros((400,0))
            bgleft_down = np.zeros((400,0))
            bgright_down = np.zeros((400,0))

            for i in range(len(result)):
                if(type(result[i]) == type(None)):
                    continue
                else:
                    blu,bru,bld,brd = result[i]
                    bgleft_up=np.column_stack((bgleft_up,blu))
                    bgright_up=np.column_stack((bgright_up,bru))
                    bgleft_down=np.column_stack((bgleft_down,bld))
                    bgright_down=np.column_stack((bgright_down,brd))

            print('Elapsed time for background estimation: ' + str(np.round((time.time() - t_background_start) / 60, 3)) + ' min')
            return bgleft_up, bgright_up, bgleft_down, bgright_down

    def getMean(self, df, mask='0'):

        if mask != '0':
            mask = mask.split(':')
            mask_chr = mask[0]
            mask_start = int(mask[1].split('-')[0])
            mask_end = int(mask[1].split('-')[1])

        def iterate_idx(i):
            np.seterr(divide='ignore', invalid='ignore')
            chr = df['chr'].iloc[i]
            x_start_index = int(df['pos1'].iloc[i])
            x_end_index = int(df['pos2'].iloc[i])
            y_start_index = int(df['pos3'].iloc[i])
            y_end_index = int(df['pos4'].iloc[i])
            region1=str(chr)+':'+str(y_start_index)+'-'+str(y_end_index)
            region2=str(chr)+':'+str(x_start_index)+'-'+str(x_end_index)
            center = self.unbalLib.fetch(region1, region2)

            centerm = np.nanmean(center)
            centersum=np.nansum(center)
            return i,centerm,centersum

        # Scoring
        ### Sobel-like operators
        nrows=len(df)
        result = Parallel(n_jobs=self.core)(delayed(iterate_idx)(i) for i in tqdm(range(nrows)))
        listM = [0 for x in range(0,nrows)]
        listS = [0 for x in range(0, nrows)]
        for r in range(len(result)):
            i,MEAN,SUM = result[r]
            listM[i] = MEAN
            listS[i] = SUM
        return listM,listS

    def pvalue(self, bgleft_up, bgright_up, bgleft_down, bgright_down, df):
        background_size = 50000 / self.resol
        background_size = int(background_size)

        np.seterr(divide='ignore', invalid='ignore')
        PVAL = []
        dfsize = len(df)
        for i in range(dfsize):
            chr = df['chr'].iloc[i]
            chr = str(chr)
            chrLen = self.chromnames2sizes[chr]
            pos1 = df['pos1'].iloc[i]
            pos2 = df['pos2'].iloc[i]
            pos3 = df['pos3'].iloc[i]
            pos4 = df['pos4'].iloc[i]
            leftmost = pos1 - background_size * self.resol
            rightmost = pos2 + background_size * self.resol
            if leftmost < 1:
                leftmost = 1
            if rightmost > chrLen:
                rightmost = chrLen
            cd1 = chr + ":" + str(leftmost) + "-" + str(rightmost)
            cd2 = chr + ":" + str(pos3) + "-" + str(pos4)
            mat = self.unbalLib.fetch(cd2, cd1)
            mat_center = mat[:, background_size:(-1*background_size)]
            mat_left = mat[:, :background_size]
            mat_right = mat[:, (-1*background_size):]

            mat_center = nantozero(mat_center)
            mat_left = nantozero(mat_left)
            mat_right = nantozero(mat_right)

            center = np.mean(mat_center, axis=1)
            left = np.mean(mat_left, axis=1)
            right = np.mean(mat_right, axis=1)

            left_diff = np.subtract(center, left)
            right_diff = np.subtract(center, right)

            pvalues = []

            x1 = (pos1 - 1) / self.resol
            x2 = pos2 / self.resol
            y1 = (pos3 - 1) / self.resol
            y2 = pos4 / self.resol

            for j in range(len(center)):
                if x1 == y1:  # downward stripe
                    difference = j
                    if difference >= 400:
                        difference = 399
                    difference = int(difference)
                    bleft = bgleft_down[difference,:]
                    bright = bgright_down[difference,:]
                elif x2 == y2:  # upward stripe
                    difference = y2 - y1 - j - 1
                    if difference >= 400:
                        difference = 399
                    difference = int(difference)
                    bleft = bgleft_up[difference,:]
                    bright = bgright_up[difference,:]

                p1 = len(np.where(bleft >= left_diff[j])[0]) / len(bleft[np.where(~np.isnan(bleft))])
                p2 = len(np.where(bright >= right_diff[j])[0]) / len(bright[np.where(~np.isnan(bright))])
                pval = max(p1, p2)
                if pval == 0:
                    pval = 1/len(bleft)
                pvalues.append(pval)
            PVAL.append(np.median(pvalues))
        return PVAL

    def pvalue_test(self, bgleft, bgright, df):
        np.seterr(divide='ignore', invalid='ignore')
        PVAL = []
        dfsize = len(df)
        for i in range(dfsize):
            chr = df['chr'].iloc[i]
            pos1 = df['pos1'].iloc[i]
            pos2 = df['pos2'].iloc[i]
            pos3 = df['pos3'].iloc[i]
            pos4 = df['pos4'].iloc[i]
            medpixel = df['medpixel'].iloc[i]
            cd1 = chr + ":" + str(pos1 - 5 * self.resol) + "-" + str(pos2 + 5 * self.resol)
            cd2 = chr + ":" + str(pos3) + "-" + str(pos4)
            mat = self.unbalLib.fetch(cd2, cd1)
            mat_center = mat[:, 5:-5]
            mat_left = mat[:, :5]
            mat_right = mat[:, -5:]

            mat_center = nantozero(mat_center)
            mat_left = nantozero(mat_left)
            mat_right = nantozero(mat_right)

            center = np.mean(mat_center, axis=1)
            left = np.mean(mat_left, axis=1)
            right = np.mean(mat_right, axis=1)

            center = signal.medfilt(center, 3)
            left = signal.medfilt(left, 3)
            right = signal.medfilt(right, 3)

            left_diff = np.subtract(center, left)
            right_diff = np.subtract(center, right)

            pcrit = []
            pvalues = []
            nomore = 0

            x1 = (pos1 - 1) / self.resol
            x2 = pos2 / self.resol
            y1 = (pos3 - 1) / self.resol
            y2 = pos4 / self.resol

            for j in range(len(center)):
                if x1 == y1:  # downward stripe
                    difference = j
                elif x2 == y2:  # upward stripe
                    difference = y2 - y1 - j
                difference = int(difference)
                bleft = bgleft[difference, :]
                bright = bgright[difference, :]
                p1 = len(np.where(bleft >= left_diff[j])[0]) / len(bleft)
                p2 = len(np.where(bright >= right_diff[j])[0]) / len(bright)
                pval = max(p1, p2)
                if pval == 0:
                    pval = 0.001
                pvalues.append(pval)
            pvalues_flip = pvalues[::-1]
            pvalue_smooth = signal.medfilt(pvalues, 3)
            pvalues_flip_smooth = signal.medfilt(pvalues_flip, 3)
            N = len(pvalue_smooth)
            for j in range(len(center)):
                if x1 == y1:
                    k = j
                elif x2 == y2:
                    k = N - j - 1
                if center[k] > medpixel and pvalue_smooth[k] < 0.5:
                    pcrit.append(1)
                elif center[k] > medpixel and pvalue_smooth[k] >= 0.5:
                    if nomore != 0:
                        pcrit.append(0)
                    else:
                        pcrit.append(1)
                elif center[k] < medpixel:
                    if nomore == 0:
                        nomore = 1
                    pcrit.append(0)

            zeros = [i for i in range(len(pcrit)) if pcrit[i] == 0]

            PVAL.append(np.median(pvalues))
        return PVAL

    def scoringstripes(self, df, expecVal, mask='0'):
        background_size = 50000/self.resol
        background_size = int(background_size)

        if mask != '0':
            mask = mask.split(':')
            mask_chr = mask[0]
            mask_start = int(mask[1].split('-')[0])
            mask_end = int(mask[1].split('-')[1])
            mask_x_start = int(mask_start / self.resol)
            mask_x_end = int(mask_end / self.resol)

        def masking(matrix ,mask_index_start, mask_index_end, start_index, end_index, direc):
            matrix = matrix.astype('float')
            relative_mask_index_start = mask_index_start - start_index
            mask_size = mask_index_end - mask_index_start
            relative_mask_index_end = relative_mask_index_start + mask_size
            L = end_index - start_index +1
            idx = [x for x in range(relative_mask_index_start, relative_mask_index_end+1) if x in range(L)]
            if len(idx) == 0:
                return matrix
            else:
                if direc == 1:
                    for x in range(matrix.shape[0]):
                        for y in idx:
                            matrix[x,y] = np.nan
                if direc == 2:
                    for x in idx:
                        for y in range(matrix.shape[1]):
                            matrix[x,y] = np.nan

            return matrix

        def expecMatrix(exval, x_start_index, x_end_index, y_start_index, y_end_index):
            x_seq = range(x_start_index, x_end_index)
            y_seq = range(y_start_index, y_end_index)
            x_seq = np.array(x_seq)
            y_seq = np.array(y_seq)

            index_matrix = [[abs(x_seq - y_seq[i])] for i in range(len(y_seq))]
            index_matrix = np.array(index_matrix)
            expec_matrix = np.empty(shape = (index_matrix.shape[0],index_matrix.shape[2]))

            for x in range(len(x_seq)):
                for y in range(len(y_seq)):
                    idx = index_matrix[y][0][x]
                    if idx >= 400:
                        idx = 399
                    expec_matrix[y,x] = exval[idx]

            del index_matrix
            return expec_matrix

        def iterate_idx(i,is_mask,exval):
            np.errstate(divide='ignore', invalid='ignore')
            xs = df['pos1'].iloc[i]
            xe = df['pos2'].iloc[i]
            ys = df['pos3'].iloc[i]
            ye = df['pos4'].iloc[i]
            #x_start_index = int((xs - 1) / self.resol)
            x_start_index = int((xs) / self.resol) # to be deleted
            x_end_index = int(xe / self.resol)

            #y_start_index = int((ys - 1) / self.resol)
            y_start_index = int((ys) / self.resol)# to be deleted
            y_end_index = int(ye / self.resol)

            leftmost = x_start_index - background_size
            rightmost = x_end_index + background_size
            if leftmost < 1:
                leftmost = 1
            if rightmost > chrom_bin_size:
                rightmost = chrom_bin_size-1

            x_coord = str(c) + ":" + str(int(xs)) + '-' + str(int(xe))
            y_coord = str(c) + ":" + str(int(ys)) + '-' + str(int(ye))
            center_obs = self.unbalLib.fetch(y_coord, x_coord)
            center_exp = expecMatrix(exval, x_start_index, x_end_index, y_start_index, y_end_index)
            center_exp += .00000001
            center = np.divide(center_obs , center_exp)

            x_coord = str(c) + ":" + str(int(leftmost * self.resol)) + '-' + str(int(x_start_index * self.resol))
            left_obs = self.unbalLib.fetch(y_coord, x_coord)
            left_exp = expecMatrix(exval, leftmost, x_start_index, y_start_index, y_end_index)
            left_exp += .00000001
            left = np.divide(left_obs, left_exp)

            x_coord = str(c) + ":" + str(int(x_end_index * self.resol)) + '-' + str(int(rightmost * self.resol))
            right_obs = self.unbalLib.fetch(y_coord,x_coord)
            right_exp = expecMatrix(exval, x_end_index + 1, rightmost +1, y_start_index, y_end_index)
            right_exp += .00000001
            right = np.divide(right_obs, right_exp)

            if is_mask and mask_start > np.min(xs-50000, ys) and mask_start < np.max(xe+50000,ye):
                center = masking(center, mask_x_start, mask_x_end, x_start_index, x_end_index, 1)
                center = masking(center, mask_x_start, mask_x_end, y_start_index, y_end_index, 2)
                left = masking(left, mask_x_start, mask_x_end, leftmost, x_start_index, 1)
                left = masking(left, mask_x_start, mask_x_end, y_start_index, y_end_index, 2)
                right = masking(right, mask_x_start, mask_x_end, x_end_index+1, rightmost, 1)
                right = masking(right, mask_x_start, mask_x_end, y_start_index, y_end_index, 2)

            centerm = np.mean(center, axis=1)
            leftm = np.mean(left, axis=1)
            rightm = np.mean(right, axis=1)

            centerTotal = np.sum(center[np.where(~np.isnan(center))[0]])
            centerMean = np.mean(center[np.where(~np.isnan(center))[0]])

            g_xl = stats.elementwise_product_sum(K_xl, leftm, centerm)
            g_xr = stats.elementwise_product_sum(K_xr, centerm, rightm)
            g_y = stats.elementwise_product_sum(K_y, leftm, centerm, rightm)
            g_x = np.minimum(g_xl, g_xr)
            diff = [a - b for a, b in zip(g_x, g_y)]
            diff = [x for x in diff if x >= 0 or x < 0]
            avgdiff = np.mean(diff)
            g = np.nanmedian(centerm) * avgdiff
            g = float(g)
            return i,g,centerMean, centerTotal

        # Scoring
        ### Sobel-like operators
        listg = [0 for x in range(df.shape[0])]
        listMean = [0 for x in range(df.shape[0])]
        listTotal = [0 for x in range(df.shape[0])]

        K_xl = np.array([[-1, 1], [-2, 2], [-1, 1]])
        K_xr = np.array([[1, -1], [2, -2], [1, -1]])
        K_y = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]])
        chrset = list(set(df['chr']))
        for c in chrset:
            is_mask = False
            if mask!='0':
                is_mask = ( mask_chr == c )
            print(c)
            idx = np.where(df['chr'] == c)[0].tolist()
            chrom_idx = self.chromnames.index(c)
            chrom_bin_size = np.ceil(self.chromsizes[chrom_idx] / self.resol)
            chrom_bin_size = int(chrom_bin_size)
            exval = expecVal[c]

            result = Parallel(n_jobs=self.core)(delayed(iterate_idx)(i,is_mask,exval) for i in tqdm(idx))
            for r in range(len(result)):
                i,g,cM,cT = result[r]
                listg[i] = g
                listMean[i] = cM
                listTotal[i] = cT
        return listg,listMean,listTotal

    def extract(self, MP, index, perc , bgleft_up, bgright_up, bgleft_down, bgright_down):
        with np.errstate(divide='ignore', invalid='ignore'):
            def search_frame(idx):
                start = idx * 200 - 100
                end = (idx + 1) * 200 + 99
                if end >= rowsize:
                    end = rowsize - 1
                if idx == 0:
                    start = 0

                framesize = end - start + 1
                start_array = [(start + j) * self.resol + 1 for j in range(framesize)]
                end_array = [s + self.resol - 1 for s in start_array]
                last = end_array[-1]
                if last >= self.chromsizes[chridx]:
                    end_array[-1] = self.chromsizes[chridx]
                locus = chr + str(":") + str(start_array[0]) + str('-') + str(end_array[-1])
                D = self.unbalLib.fetch(locus, locus)
                D = nantozero(D)

                # Remove rows and columns containing only zero.
                colsum = np.sum(D, axis=0)
                # rowsum = np.sum(D, axis=1) # rowsum == colsum
                nonzero_idx = np.where(colsum != 0)  # data type: tuple
                nonzero_idx = nonzero_idx[0]
                framesize = len(nonzero_idx)
                if framesize > 10:
                    D = D[np.ix_(nonzero_idx, nonzero_idx)]
                    start_array = [start_array[s] for s in nonzero_idx]
                    end_array = [end_array[s] for s in nonzero_idx]
                    temp_res = self.StripeSearch(D, idx, start, end, M, perc, chr, framesize, start_array, end_array)
                    return(temp_res)
                else:
                    return 0

            result = pd.DataFrame(columns=['chr', 'pos1', 'pos2', 'chr2', 'pos3', 'pos4', 'length', 'width', 'total', 'Mean',
                         'maxpixel', 'num', 'start', 'end', 'x', 'y', 'h', 'w', 'medpixel'])
            for chridx in range(len(self.chromnames)):
                t_chr_search = time.time()
                chr = self.chromnames[chridx]
                print('Chromosome: ' + str(chr) + " / Maximum pixel: " + str(round(perc*100,3))+"%")
                #cfm = self.unbalLib.fetch(chr)  # cfm : contact frequency matrix
                #cfm = np.round(cfm,1)
                #M = np.quantile(a=cfm[cfm > 0], q=mp, interpolation='linear')
                M = MP[chr][index]
                rowsize = int(np.ceil(self.chromsizes[chridx] / self.resol))
                #rowsize = cfm.shape[0]
                #del cfm

                nframes = math.ceil(rowsize / 200)
                results = Parallel(n_jobs = self.core)(delayed(search_frame)(i) for i in tqdm(range(nframes)))

                for n in range(len(results)):
                    if type(results[n]) == int:
                        pass
                    else:
                        result = result.append(results[n])

            res = self.RemoveRedundant(result, 'size')
            res = res.reset_index(drop=True)

            # Stripe filtering and scoring
            # res2 = self.scoringstripes(res)
            p = self.pvalue(bgleft_up, bgright_up, bgleft_down, bgright_down, res)
            #s = self.scoringstripes(res)
            res = res.assign(pvalue=pd.Series(p))
            #res = res.assign(Stripiness=pd.Series(s[5]))

        return res

    def StripeSearch(self, submat, num, start, end, M, perc, chr, framesize, start_array, end_array):
        with np.errstate(divide='ignore', invalid='ignore'):
            res_chr = [];
            res_pos1 = [];
            res_pos2 = [];
            res_pos3 = [];
            res_pos4 = [];
            res_length = [];
            res_width = [];
            res_total = [];
            res_Mean = [];
            res_maxpixel = [];
            res_num = [];
            res_start = [];
            res_end = [];
            res_x = [];
            res_y = [];
            res_h = [];
            res_w = [];
            res_medpixel = []
            mp = perc
            medpixel = np.quantile(a=submat[submat > 0], q=0.5)
            st = start
            en = end
            S = framesize
            red = np.ones((framesize, framesize)) * 255
            blue = 255 * (M - submat) / M
            blue[np.where(blue < 0)] = 0
            green = blue

            img = cv.merge((red / 255, green / 255, blue / 255))
            img = np.clip(img, a_min=0, a_max=1)
            # plt.subplot(111),plt.imshow(img),plt.title('original'), plt.show()

            for b in np.arange(0.5, 1.01, 0.1):  # b: brightness parameter
                test_column = []
                end_points = []
                start_points = []
                updown = []  # up = 1, down = 2
                adj = ImageProcessing.imBrightness3D(img, In=([0.0, 0.0, 0.0], [1.0, b, b]),
                                                     Out=([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))
                # plt.subplot(111), plt.imshow(adj), plt.title('Brightened'), plt.show()
                kernel = np.ones((3, 3)) / 9
                blur = cv.filter2D(adj, -1, kernel)
                blur = np.clip(blur, a_min=0, a_max=1)
                # plt.subplot(111), plt.imshow(blur), plt.title('Blurry'), plt.show()
                gray = cv.cvtColor(np.float32(blur), cv.COLOR_RGB2GRAY)
                # gray = np.uint8(gray*255)
                # edges = ImageProcessing.Canny(gray, threshold = 0.3)
                edges = feature.canny(gray, sigma=self.canny)
                # plt.subplot(111), plt.imshow(edges, cmap='gray'), plt.title('Canny edge detection'), plt.show()
                vert = ImageProcessing.verticalLine(edges, L=60, H=120)
                # plt.subplot(111), plt.imshow(vert, cmap='gray'), plt.title('Vertical line detection'), plt.show()
                LL = []
                for c in range(S):
                    # t1 = time.time()
                    line_length, END = ImageProcessing.block(vert, c)
                    # print(time.time()-t1)
                    LL.append(line_length)
                    above = min(c, END)
                    bottom = max(c, END)
                    seq = list(range(above, bottom + 1, 1))

                    if line_length > self.minH and sum(vert[seq, c]) != 0:
                        test_column.append(c)
                        end_points.append(END)
                        start_points.append(c)
                        if END > c:
                            updown.append(2)
                        else:
                            updown.append(1)

                Pair = []
                MIN_vec = []
                MAX_vec = []
                for ud in [1, 2]:
                    testmat = np.zeros((S, S), dtype=np.uint8)
                    udidx = [i for i in range(len(updown)) if updown[i] == ud]
                    for c in udidx:
                        st = test_column[c]
                        en = end_points[c]
                        if ud == 1:
                            en_temp = st
                            st = en
                            en = en_temp
                        testmat[st:en, test_column[c]] = 1
                    # line refinement
                    for r in range(S):
                        vec = testmat[r, :]
                        K1 = vec[1:S] > vec[0:(S - 1)]
                        K2 = vec[1:S] < vec[0:(S - 1)]
                        st = [i + 1 for i in range(len(K1)) if K1[i]]
                        en = [i for i in range(len(K2)) if K2[i]]
                        if vec[0] == 1:
                            st.insert(0, 0)
                        if vec[S - 1] == 1:
                            en.insert(len(en), S - 1)

                        nLines = len(st)

                        for L in range(nLines):
                            origLine = edges[r, list(range(st[L], en[L] + 1, 1))]
                            SUM = sum(origLine)
                            if SUM > 0:
                                testmat[r, st[L]:en[L]] = vert[r, st[L]:en[L]]
                            else:
                                MED = int(np.round(stat.median([st[L] + en[L]]) / 2))
                                testmat[r, st[L]:en[L]] = 0
                                testmat[r, MED] = 1
                    fused_image = np.dstack((edges, testmat, testmat))

                    [_, Y] = np.where(testmat == 1)
                    uniqueCols = list(set(Y))
                    ps = pd.Series(Y)
                    counts = ps.value_counts().sort_index()
                    counts = counts.to_frame(name='length')
                    start_points_ud = [start_points[i] for i in udidx]
                    end_points_ud = [end_points[i] for i in udidx]
                    intersectidx = [i for i in range(len(start_points_ud)) if start_points_ud[i] in uniqueCols]
                    start_points_ud = [start_points_ud[i] for i in intersectidx]
                    end_points_ud = [end_points_ud[i] for i in intersectidx]

                    counts['end_points'] = end_points_ud

                    counts = counts[counts['length'] >= 3]
                    nrow = counts.shape[0]
                    meanX = []
                    Continuous = []
                    isContinue = False
                    Len = []

                    for c in range(nrow - 1):
                        Current = counts.index[c]
                        Next = counts.index[c + 1]

                        if Next - Current == 1 and isContinue:
                            Continuous.append(Next)
                            Len.append(counts.iloc[c + 1]['length'])
                        elif Next - Current == 1 and not isContinue:
                            Continuous = [Current, Next]
                            Len = [counts.iloc[c]['length'], counts.iloc[c + 1]['length']]
                            isContinue = True
                        elif Next - Current != 1 and not isContinue:
                            Continuous = [Current]
                            Len = [counts.iloc[c]['length']]
                            Len = [a / sum(Len) for a in Len]
                            isContinue = False
                            temp = sum([a * b for a, b in zip(Continuous, Len)])
                            meanX.append(np.round(temp))
                        else:
                            Len = [a / sum(Len) for a in Len]
                            temp = sum([a * b for a, b in zip(Continuous, Len)])
                            meanX.append(np.round(temp))
                            Continuous = [Current]
                            Len = [counts.iloc[c]['length']]
                            isContinue = False
                    Len = [a / sum(Len) for a in Len]
                    temp = sum([a * b for a, b in zip(Continuous, Len)])
                    meanX.append(np.round(temp))

                    X = list(set(meanX))
                    X.sort()
                    Xsize = len(X)

                    for c in range(Xsize - 1):
                        n = int(X[c])
                        m = int(X[c + 1])
                        st1 = np.where(testmat[:, n] == 1)[0]
                        en1 = st1.max()
                        st1 = st1.min()
                        st2 = np.where(testmat[:, m] == 1)[0]
                        en2 = st2.max()
                        st2 = st2.min()
                        '''
                        max_width_dist1.append(abs(m - n))
                        if c == 0:
                            max_width_dist2.append(abs(m - n))
                        else:
                            l = int(X[c - 1])
                            minw = min(abs(m - n), abs(n - l))
                            if max_width_dist2[-1] == minw:
                                continue
                            else:
                                max_width_dist2.append(minw)
                        '''
                        if abs(m - n) > 1 and abs(m - n) <= self.maxW:
                            Pair.append((n, m))
                            [a1, _] = np.where(testmat[:, range(max(0, n - 1), min(n + 2, S), 1)] == 1)
                            MIN1 = a1.min()
                            MAX1 = a1.max()
                            [a2, _] = np.where(testmat[:, range(max(0, m - 1), min(m + 2, S), 1)] == 1)
                            MIN2 = a2.min()
                            MAX2 = a2.max()

                            MIN = min(MIN1, MIN2)
                            MAX = max(MAX1, MAX2)

                            if ud == 1:
                                MAX = X[c + 1]
                            else:
                                MIN = X[c]
                            MIN_vec.append(MIN)
                            MAX_vec.append(MAX)
                PairSize = len(Pair)

                for c in range(PairSize):
                    x = Pair[c][0]
                    y = int(MIN_vec[c])
                    w = int(Pair[c][1] - Pair[c][0] + 1)
                    h = int(MAX_vec[c] - MIN_vec[c] + 1)

                    res_chr.append(chr)
                    res_pos1.append(start_array[x])
                    res_pos2.append(end_array[x + w - 1])
                    res_pos3.append(start_array[y])
                    res_pos4.append(end_array[y + h - 1])
                    res_length.append(end_array[y + h - 1] - start_array[y] + 1)
                    res_width.append(end_array[x + w - 1] - start_array[x] + 1)
                    res_total.append(submat[y:(y + h), x:(x + w)].sum())
                    res_Mean.append(submat[y:(y + h), x:(x + w)].sum() / h / w)
                    res_maxpixel.append(str(mp * 100) + '%')
                    res_num.append(num)
                    res_start.append(start)
                    res_end.append(end)
                    res_x.append(x)
                    res_y.append(y)
                    res_h.append(h)
                    res_w.append(w)
                    res_medpixel.append(medpixel)

            result = pd.DataFrame(
                {'chr': res_chr, 'pos1': res_pos1, 'pos2': res_pos2, 'chr2': res_chr, 'pos3': res_pos3, 'pos4': res_pos4,
                 'length': res_length, 'width': res_width, 'total': res_total, 'Mean': res_Mean,
                 'maxpixel': res_maxpixel, 'num': res_num, 'start': res_start, 'end': res_end,
                 'x': res_x, 'y': res_y, 'h': res_h, 'w': res_w, 'medpixel': res_medpixel})

            result = self.RemoveRedundant(result, 'size')

        return result

    def RemoveRedundant(self, df, by):
        def clean(n):
            delidx=[]
            n_idx = np.where(subdf['num'] == n)[0]
            n2_idx = np.where(subdf['num'] == n + 1)[0]
            n_idx = np.concatenate((n_idx, n2_idx))
            n_idx.sort()
            L = len(n_idx)
            for i in range(L - 1):
                for j in range(i + 1, L):
                    ii = c_idx[n_idx][i]
                    jj = c_idx[n_idx][j]

                    A_x_start = list_pos1[ii]
                    A_x_end = list_pos2[ii]
                    A_y_start = list_pos3[ii]
                    A_y_end = list_pos4[ii]

                    B_x_start = list_pos1[jj]
                    B_x_end = list_pos2[jj]
                    B_y_start = list_pos3[jj]
                    B_y_end = list_pos4[jj]

                    int_x = range(max(A_x_start, B_x_start), min(A_x_end, B_x_end) + 1)
                    int_y = range(max(A_y_start, B_y_start), min(A_y_end, B_y_end) + 1)

                    s_x = len(int_x) / min(A_x_end - A_x_start, B_x_end - B_x_start)
                    s_y = len(int_y) / min(A_y_end - A_y_start, B_y_end - B_y_start)

                    if s_x > 0.2 and s_y > 0.2:
                        if by == 'size':
                            if list_h[ii] * list_w[ii] <= list_h[jj] * list_w[jj]:
                                delidx.append(ii)
                            else:
                                delidx.append(jj)
                        elif by == 'score':
                            if list_stri[ii] <= list_stri[jj]:
                                delidx.append(ii)
                            else:
                                delidx.append(jj)
                        else:
                            if list_pval[ii] > list_pval[jj]:
                                delidx.append(ii)
                            else:
                                delidx.append(jj)
            return delidx

        if by != 'size' and by != 'score' and by != 'pvalue':
            raise ValueError('"by" should be one of "size", "pvalue" and "score"')
        df_size = df.shape
        row_size = df_size[0]
        if row_size == 0:
            return df
        else:
            delobj = [True for i in range(row_size)]
            list_chr = df['chr']
            list_pos1 = df['pos1'].tolist()
            list_pos2 = df['pos2'].tolist()
            list_pos3 = df['pos3'].tolist()
            list_pos4 = df['pos4'].tolist()
            list_h = df['h'].tolist()
            list_w = df['w'].tolist()
            if by == 'score':
                list_stri = df['Stripiness'].tolist()
            if by == 'pvalue':
                list_pval = df['pvalue'].tolist()
            unique_chr = list(set(list_chr))

            for c in unique_chr:
                c_idx = np.where(list_chr == c)[0]
                subdf = df.iloc[c_idx]
                unique_num = list(set(subdf['num']))
                unique_num.sort()
                res = Parallel(n_jobs=self.core)(delayed(clean)(n) for n in unique_num)
                for i in range(len(res)):
                    for j in res[i]:
                        delobj[j] = False

        idx = [a for a in range(row_size) if delobj[a]]
        result = df.iloc[idx]
        return result

    def selectColumn(self, df):
        list_chr = df['chr']
        list_pos1 = df['pos1'].tolist()
        list_pos2 = df['pos2'].tolist()
        list_pos3 = df['pos3'].tolist()
        list_pos4 = df['pos4'].tolist()

        if str(list_chr[0])[0] != 'c':
            if self.chromnames[0][0] == 'c':
                list_chr = ['chr'+str(x) for x in list_chr]
        MAX_POS=[]
        nrow = df.shape[0]
        for i in range(nrow):
            chr = list_chr[i]
            pos1 = list_pos1[i]
            pos2 = list_pos2[i]
            pos3 = list_pos3[i]
            pos4 = list_pos4[i]

            x_str = str(chr)+':'+str(pos1)+'-'+str(pos2)
            y_str = str(chr)+':'+str(pos3)+'-'+str(pos4)
            mat = self.unbalLib.fetch(x_str,y_str)
            average = np.mean(mat,axis=1)
            which_max = np.argmax(average)
            max_pos = pos1 + self.resol * which_max
            MAX_POS.append(max_pos)


        return MAX_POS


def nantozero(nparray):
    where_are_nans = np.isnan(nparray)
    nparray[where_are_nans] = 0
    return nparray


def extract_from_matrix(matrix, x_start, x_end, y_start, y_end, mask_x_start = 0, mask_x_end = 0):
    x = list(range(x_start, x_end))
    x_mask = list(range(mask_x_start, mask_x_end))
    x = [i for i in x if i not in x_mask]
    xlen = len(x)
    y = list(range(y_start, y_end))
    y = [i for i in y if i not in x_mask]
    ylen = len(y)
    result = np.empty((ylen, xlen), dtype=float)
    for i in range(ylen):
        for j in range(xlen):
            result[i][j] = matrix[y[i]][x[j]]


    return result


