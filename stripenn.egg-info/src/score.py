import pandas as pd
import cooler
import numpy as np
from src import getStripe

def getScore(cool, coordinates,norm,numcores,out,mask='0'):
    # Just open without header
    table = pd.read_csv(coordinates,header=None,sep='\t')
    # test
    el = table.iloc[0]
    el1 = el[1]
    el2 = el[2]

    if type(el1) != str and type(el2) != str :
        pass
    else:
        table = pd.read_csv(coordinates,header=0, sep='\t')

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
    all_chromnames = [x for x in all_chromnames if x != "Y"]
    #all_chromsizes = Lib.chromsizes
    all_chromsizes = [Lib.chromsizes[i] for i in range(len(Lib.chromsizes)) if Lib.chromnames[i]!="Y"]
    all_chromsizes = np.array(all_chromsizes)
    chrom_remain_idx = np.where(all_chromsizes > 1000000)[0]
    all_chromnames = [all_chromnames[i] for i in chrom_remain_idx]
    all_chromsizes = all_chromsizes[chrom_remain_idx]
    chromnames = all_chromnames
    chromsizes = all_chromsizes
    unbalLib = Lib.matrix(balance=norm)
    resol = Lib._info['bin-size']
    obj = getStripe.getStripe(unbalLib, resol, 10, 8, 2.5, all_chromnames, chromnames, all_chromsizes, chromsizes,core)
    EV = getStripe.getStripe.mpmean(obj)
    bgleft_up, bgright_up, bgleft_down, bgright_down = getStripe.getStripe.nulldist(obj)
    pval = getStripe.getStripe.pvalue(obj, bgleft_up, bgright_up, bgleft_down, bgright_down, table)
    table.insert(table.shape[1],'pvalue_added',pval,True)
    MEAN, SUM = getStripe.getStripe.getMean(obj,table)
    s,MEANOE,TOTALOE = obj.scoringstripes(table, EV, mask)
    table.insert(table.shape[1],'Stripiness_added',s,True)
    table.insert(table.shape[1], "O_Mean_added", MEAN, True)
    table.insert(table.shape[1], "O_Sum_added", SUM, True)
    table.insert(table.shape[1], 'O/E_Mean_added',MEANOE,True)
    table.insert(table.shape[1], 'O/E_Total_added',TOTALOE,True)
    table.to_csv(out,sep="\t",header=True,index=False)
