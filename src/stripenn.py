import argparse
import cooler
import multiprocessing
from src import getStripe
import os
import shutil
import errno
import pandas as pd
import numpy as np
import warnings
import time
import sys

def makeOutDir(outdir):
    last = outdir[-1]
    if last != '/':
        outdir += '/'
    if os.path.exists(outdir):

        print('\n%s exists. Do you want to remove all files and save new results in this folder? [Y/n]' % (outdir))
        userinput = input()
        if userinput == 'Y' or userinput == 'y':
            print('All directories and files in %s will be deleted.' % (outdir))
            for filename in os.listdir(outdir):
                file_path = os.path.join(outdir, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('Failed to delete %s with the reason: %s' % (file_path, e))
        elif userinput == 'n' or userinput == 'N':
            print('Input another output directory. Exit.')
            quit()
        else:
            print('Type Y or n.\nExit.')
            quit()
    else:
        try:
            os.makedirs(outdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


def addlog(cool, out, norm, chrom, canny, minL, maxW, maxpixel, numcores, pvalue, mask):
    if out[-1] != '/':
        out += '/'
    outfile=open(out + "stripenn.log",'w')
    outfile.write('cool: ' + cool + '\n')
    outfile.write('out: ' + out + '\n')
    outfile.write('norm: ' + norm + '\n')
    outfile.write('chrom: ' + str(chrom) + '\n')
    outfile.write('canny: ' + str(canny) + '\n')
    outfile.write('minL: ' + str(minL) + '\n')            
    outfile.write('maxW: ' + str(maxW) + '\n')                
    outfile.write('maxpixel: ' + str(maxpixel) + '\n')
    outfile.write('num_cores: ' + str(numcores) + '\n')
    outfile.write('pvalue: ' + str(pvalue) + '\n')
    outfile.write('mask: ' + str(mask) + '\n')
    outfile.close()
    
    
def compute(cool, out, norm, chrom, canny, minL, maxW, maxpixel, numcores, pvalue, mask, slow):
    np.seterr(divide='ignore', invalid='ignore')
    t_start = time.time()
    if out[-1] != '/':
        out += '/'
    makeOutDir(out)
    addlog(cool, out, norm, chrom, canny, minL, maxW, maxpixel, numcores, pvalue, mask)
    maxpixel = maxpixel.split(',')
    maxpixel = list(map(float, maxpixel))
    chroms = chrom.split(',')
    minH = minL
    core = numcores
    pcut = pvalue
    #cool, out, norm, chroms, canny, minH, maxW, maxpixel, core, pcut = argumentParser()
    print('Result will be stored in %s' % (out))

    Lib = cooler.Cooler(cool)
    PossibleNorm = Lib.bins().columns
    if norm == 'None':
        norm = False
    elif norm == 'weight':
        norm = True
    elif norm not in PossibleNorm:
        print('Possible normalization methods are:')
        print('None')
        for n in range(3,len(PossibleNorm)):
            print(PossibleNorm[n])
        print("Invalid normalization method. Normalization method is forced to None")
        norm = False

    all_chromnames = Lib.chromnames
    all_chromsizes = Lib.chromsizes
    chrom_remain_idx = [i for i in range(len(all_chromnames)) if 'JH5' not in all_chromnames[i] and 'GL4' not in all_chromnames[i] and all_chromnames[i] != 'M' and all_chromnames[i] != "chrM" and all_chromnames[i] != "Y" and all_chromnames[i] != "chrY"]
    all_chromnames = [all_chromnames[i] for i in chrom_remain_idx]
    all_chromsizes = all_chromsizes[chrom_remain_idx]
    chromnames = all_chromnames
    chromsizes = all_chromsizes

    if len(all_chromnames) == 0:
        sys.exit("Exit: All chromosomes are shorter than 50kb.")
    warnflag = False
    if chroms[0] != 'all':
        idx = []
        for item in chroms:
            if item in all_chromnames:
                idx.append(all_chromnames.index(item))
            else:
                warnings.warn('\nThere is no chromosomes called '+str(item)+' in the provided .cool file or it is shorter than 50kb.')
                warnflag = True
        if warnflag:
            warnings.warn('\nThe possible chromosomes are: '+ ', '.join(all_chromnames))
        chromnames = chroms
        chromsizes = all_chromsizes[idx]

    unbalLib = Lib.matrix(balance=norm)
    resol = Lib.binsize
    obj = getStripe.getStripe(unbalLib, resol, minH, maxW, canny, all_chromnames, chromnames, all_chromsizes, chromsizes,core)
    print('1. Maximum pixel value calculation ...')
    if slow:
        print("1.1 Slowly estimating Maximum pixel values...")
        MP = getStripe.getStripe.getQuantile_slow(obj,Lib,chromnames,maxpixel)
    else:
        MP = getStripe.getStripe.getQuantile_original(obj,Lib,chromnames,maxpixel)
    print('2. Expected value calculation ...')
    EV = getStripe.getStripe.mpmean(obj)
    print('3. Background distribution estimation ...')
    bgleft_up, bgright_up, bgleft_down, bgright_down = getStripe.getStripe.nulldist(obj)
    print('4. Finding candidate stripes from each chromosome ...')
    result_table = pd.DataFrame(columns=['chr', 'pos1', 'pos2', 'chr2', 'pos3', 'pos4', 'length', 'width', 'total', 'Mean',
                                   'maxpixel', 'num', 'start', 'end', 'x', 'y', 'h', 'w', 'medpixel','pvalue'])
    for i in range(len(maxpixel)):
        perc = maxpixel[i]
        result = obj.extract(MP, i, perc, bgleft_up, bgright_up, bgleft_down, bgright_down)
        result_table = result_table.append(result)
#    for mp in maxpixel:
#        result = obj.extract(mp, bgleft, bgright)
#        result_table = result_table.append(result)

    result_table = getStripe.getStripe.RemoveRedundant(obj,df=result_table,by='pvalue')

    print('5. Stripiness calculation ...')
    s = obj.scoringstripes(result_table, EV, mask)
    s = s[0]
    result_table = result_table.drop(columns=['total', 'num', 'start', 'end', 'x', 'y', 'h', 'w', 'medpixel'])
    result_table.insert(result_table.shape[1],'Stripiness',s,True)
    res_filter = result_table[result_table['pvalue'] < pcut]
    res_filter = res_filter.sort_values(by=['Stripiness'], ascending=False)


    res1 = out + 'result_unfiltered.tsv'
    res2 = out + 'result_filtered.tsv'

    result_table.to_csv(res1,sep="\t",header=True,index=False)
    res_filter.to_csv(res2,sep="\t",header=True,index=False)

    print('\n'+str(round((time.time()-t_start)/60,3))+'min taken.')
    print('Check the result stored in %s' % (out))
    return 0

#if __name__ == "__main__":
#    main()

