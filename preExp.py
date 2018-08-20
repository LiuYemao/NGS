#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

def parseExp(s):
    f = "/home/ymliu/bna/11.exp/{}.exp".format(s)
    df = pd.read_table(f, sep = "\t", index_col=0)
    df.columns = ["rep1", 'rep2', 'rep3', 'm', 'g_m']
    return df
    # dic = {}
    # with open(f) as rf:
    #     rf.readline()
    #     for line in rf:
    #         ll = line.strip().split()
    #         if len(ll) == 6:
    #             dic[ll[0]] = float([-1])

    # return dic


def main(wfname):
    TD013 = parseExp("TD013");TH363 = parseExp("TH363")
    TH429 = parseExp("TH429");TH430 = parseExp("TH430")
    TH433 = parseExp("TH433");TH008 = parseExp("TH008")
    bol = parseExp("bol2bol");bra = parseExp("bra2bra")
    synfile = "/home/ymliu/ref/bna/homo_genes.txt"
    brgenes = []; bogenes = [];bnaA = [];bnaC = []
    with open(synfile) as rf:
        for line in rf:
            if '.' not in line:
                ll = line.strip().split()
                if len(ll) == 4:
                    brgenes.append(ll[0])
                    bogenes.append(ll[1])
                    bnaA.append(ll[2])
                    bnaC.append(ll[3])
    br_value = [bra.loc[i, 'g_m'] if i in bra.index else np.nan for i in brgenes]
    bo_value = [bol.loc[i, 'g_m'] if i in bol.index else np.nan for i in bogenes]
    TH008_Avalue = [TH008.loc[i, 'g_m'] if i in TH008.index else np.nan for i in bnaA]
    TH008_Cvalue = [TH008.loc[i, 'g_m'] if i in TH008.index else np.nan for i in bnaC]
    TD013_Avalue = [TD013.loc[i, 'g_m'] if i in TD013.index else np.nan for i in bnaA]
    TD013_Cvalue = [TD013.loc[i, 'g_m'] if i in TD013.index else np.nan for i in bnaC]
    TH363_Avalue = [TH363.loc[i, 'g_m'] if i in TH363.index else np.nan for i in bnaA]
    TH363_Cvalue = [TH363.loc[i, 'g_m'] if i in TH363.index else np.nan for i in bnaC]
    TH429_Avalue = [TH429.loc[i, 'g_m'] if i in TH429.index else np.nan for i in bnaA]
    TH429_Cvalue = [TH429.loc[i, 'g_m'] if i in TH429.index else np.nan for i in bnaC]
    TH430_Avalue = [TH430.loc[i, 'g_m'] if i in TH430.index else np.nan for i in bnaA]
    TH430_Cvalue = [TH430.loc[i, 'g_m'] if i in TH430.index else np.nan for i in bnaC]
    TH433_Avalue = [TH433.loc[i, 'g_m'] if i in TH433.index else np.nan for i in bnaA]
    TH433_Cvalue = [TH433.loc[i, 'g_m'] if i in TH433.index else np.nan for i in bnaC]

    dat = {"genes_Br":brgenes, "FPKM_Br":br_value, 'genes_Bo':bogenes, 'FPKM_Bo':bo_value,
           'genes_BnA':bnaA, 'FPKM_TD013_A':TD013_Avalue, 'FPKM_TH008_A':TH008_Avalue,
           'FPKM_TH363_A': TH363_Avalue,'FPKM_TH429_A':TH429_Avalue,'FPKM_TH430_A':TH430_Avalue,
           'FPKM_TH433_A': TH433_Avalue,
           'genes_BnC': bnaC, 'FPKM_TD013_C': TD013_Cvalue, 'FPKM_TH008_C': TH008_Cvalue,
           'FPKM_TH363_C': TH363_Cvalue, 'FPKM_TH429_C': TH429_Cvalue, 'FPKM_TH430_C': TH430_Cvalue,
           'FPKM_TH433_C': TH433_Cvalue}
    df = pd.DataFrame(dat)
    mydf = df.dropna(axis=0, how='any')
    df.to_csv('{}.csv'.format(wfname), index=False)
    my.to_csv('{}.allrows.csv'.format(wfname), index = False)

if __name__ == '__main__':

    main(sys.argv[1])
