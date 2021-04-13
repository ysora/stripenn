import pandas as pd
import cooler
import numpy as np
from . import getStripe

def getScore(cool, coordinates,norm,header,numcores,out):
    print(header)
    if header:
        header = 0
    else:
        header = None
    table = pd.read_csv(coordinates,header=header,sep='\t')
    table.columns=['chr','pos1','pos2','chr2','pos3','pos4']+table.columns[6:].tolist()

    core = numcores
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

    all_chromnames = Lib.chromnames
    all_chromsizes = Lib.chromsizes
    chrom_remain_idx = np.where(all_chromsizes > 500000)[0]
    all_chromnames = [all_chromnames[i] for i in chrom_remain_idx]
    all_chromsizes = all_chromsizes[chrom_remain_idx]
    chromnames = all_chromnames
    chromsizes = all_chromsizes
    unbalLib = Lib.matrix(balance=norm)
    resol = Lib._info['bin-size']
    obj = getStripe.getStripe(unbalLib, resol, 10, 8, 2.5, all_chromnames, chromnames, all_chromsizes, chromsizes,core)
    bgleft, bgright = getStripe.getStripe.nulldist(obj)
    pval = getStripe.getStripe.pvalue(obj, bgleft, bgright, table)
    table.insert(table.shape[1],'pvalue_added',pval,True)

    s = obj.scoringstripes(table)
    table.insert(table.shape[1],'Stripiness_added',s,True)
    table.to_csv(out,sep="\t",header=True,index=False)
