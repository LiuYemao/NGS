#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
========
使用方法：
1. 默认阈值筛选SNP（maf：不小于0.05； het：不大于0.05；miss，保留全部）：
python processVCF.py -i <读入的VCF文件名> -o <输出VCF文件名>

2. 保留het（大于0.06，小于0.08），maf（大于0.05 小于0.07）的SNP，miss大于0.15的SNP：
python processVCF.py -i <读入的VCF文件名> -o <输出VCF文件名> -e 0.06 -t 0.08 -a 0.05 -x 0.07 -m 0.15 

3. 提取出部分样本的基因型
首先准备好一个样本名的文件，一行一个样本名
python processVCF.py -i <读入的VCF文件名> -o <输出VCF文件名> -I <含有样本的文件>

========

@author: ymliu
"""

from __future__ import division
import sys
import numpy as np

chrAs = ["chrA01","chrA02","chrA03","chrA04","chrA05","chrA06","chrA07","chrA08","chrA09","chrA10"]
chrCs = ["chrC01","chrC02","chrC03","chrC04","chrC05","chrC06","chrC07","chrC08","chrC09"]
chrs = chrAs + chrCs
chrs_2 = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

def listFreq(L):
    dic = {}
    for l in set(L):
        dic[l] = L.count(l)
    return dic

def openGzipFile(f):
    import gzip
    return gzip.open(f)

def checkFileType(line):
    if "/" in line and "|" in line:
        print "The VCF file is not effective. It should be a phased or unphased file, not a mixed one." 
        return False
    if "/" in line and "|" not in line:
        return "/"
    elif "|" in line and "/" not in line:
        return "|"
    else:
        print "The site in genotype VCF file must be biallelic site."
        return False

def preData(readInFile, minHetThreshold, maxHetThreshold, minMAFThreshold, maxMAFThreshold,\
            maxMissThreshold, writeOutFile = None, writeHetFile = None, filterIndividualFile = None,\
            writeMafFile = None, writeReportFile = None, writeDPFile = None, convert=False):     
    
    minHetThreshold = float(minHetThreshold)
    maxHetThreshold = float(maxHetThreshold)
    minMAFThreshold = float(minMAFThreshold)
    maxMAFThreshold = float(maxMAFThreshold)
    maxMissThreshold = float(maxMissThreshold)
          
    if writeOutFile != None: writeOut = open(writeOutFile, 'w')
    if writeHetFile != None: 
        writeHet = open(writeHetFile, 'w')
        writeHet.write('\t'.join(['snp','heterogeneity'])+'\n')
    if writeMafFile != None: 
        writeMaf = open(writeMafFile, 'w')
        writeMaf.write('\t'.join(['snp','maf'])+'\n')
    if writeDPFile != None: 
        writeDP = open(writeDPFile, 'w')
        writeDP.write('\t'.join(['snp','depth'])+'\n')
    if writeReportFile != None:
        writeReport = open(writeReportFile, "w")
        writeReport.write('\t'.join(['snp', 'maf', 'het', 'miss', 'DP'])+'\n')
    if filterIndividualFile != None:
        with open(filterIndividualFile) as indFile:
            samplesDic = {}
            indivs = []
            for line in indFile:
                ll = line.rstrip('\r\n').split()
                if len(ll)==2:
                    samplesDic[ll[0]] = ll[1]
                    indivs.append(ll[0])
                elif len(ll) == 1:
                    samplesDic[ll[0]] = ll[0]
                    indivs.append(ll[0])
    if readInFile.endswith(".vcf.gz"):
        readIn = openGzipFile(readInFile)
    elif readInFile.endswith(".vcf"):
        readIn = open(readInFile)
    else:
        sys.exit()        
    for line in readIn:           
        if line.upper().startswith("#CHROM"):
            if filterIndividualFile != None:
                tmplist = line.rstrip('\r\n').split()
                indexes = [tmplist.index(i) for i in indivs]                
                if writeOutFile != None:
                    writeOut.write('\t'.join(tmplist[:9]\
                        +[samplesDic[tmplist[i]] for i in indexes])+'\n')
            else:
                if writeOutFile != None:
                    writeOut.write(line)
        elif line.startswith("#"):
            if writeOutFile != None:
                writeOut.write(line)
        elif not line.startswith('#'):
            
            if filterIndividualFile != None:
                ll = []
                tmplist = line.rstrip('\r\n').split()
                ll = tmplist[:9]
                for j in indexes:
                    ll.append(tmplist[j])
                if writeDPFile != None or writeReportFile != None:
                    DPs = []
                    try:
                        DP_index = ll[8].split(':').index("DP")
                    except:
                        sys.exit("The vcf file doesn't contain any DP information! \
                            Please check the file.")                    
                for i,l in enumerate(ll):
                    if ':' in l:
                        tmplistl = l.split(":")
                        if writeDPFile != None or writeReportFile != None:
                            if len(tmplistl) > DP_index and i>=9 and \
                            tmplistl[DP_index] != ".":
                                DPs.append(int(tmplistl[DP_index]))
                        ll[i] = tmplistl[0]                   
            else:
                ll = line.rstrip('\r\n').split()
                if writeDPFile != None or writeReportFile != None:
                    DPs = []
                    try:
                        DP_index = ll[8].split(':').index("DP")
                    except:
                        sys.exit("The vcf file doesn't contain any DP information! \
                            Please check the file.")
                for i,l in enumerate(ll):
                    if ':' in l:
                        tmplistl = l.split(":")
                        if writeDPFile != None or writeReportFile != None:
                            if len(tmplistl) > DP_index and i>=9 and \
                            tmplistl[DP_index] != ".":
                                DPs.append(int(tmplistl[DP_index]))
                        ll[i] = tmplistl[0]

            if convert != False:
                chrDic = {'chrA01':'1', 'chrA02':'2', 'chrA03':'3', 'chrA04':'4', 'chrA05':'5',\
                          'chrA06':'6', 'chrA07':'7', 'chrA08':'8', 'chrA09':'9', 'chrA10':'10',\
                          'chrC01':'11', 'chrC02':'12', 'chrC03':'13', 'chrC04':'14', 'chrC05':'15',
                          'chrC06':'16', 'chrC07':'17', 'chrC08':'18', 'chrC09':'19',\
                          'chrA01_random':"20", 'chrA02_random':"21", 'chrA03_random':"22",\
                          'chrA04_random':"23", 'chrA05_random':"24", 'chrA06_random':"25",\
                          'chrA07_random':"26", 'chrA08_random':"27", 'chrA09_random':"28",\
                          'chrA10_random':"29", 'chrAnn_random':"30", 'chrC01_random':"31",\
                          'chrC02_random':"32", 'chrC03_random':"33", 'chrC04_random':"34",\
                          'chrC05_random':"35", 'chrC06_random':"36", 'chrC07_random':"37",\
                          'chrC08_random':"38", 'chrC09_random':"39", 'chrCnn_random':"40",\
                          'chrUnn_random':"41"}
                ll[0] = chrDic[ll[0]]                  
            
            if len(ll[3])+len(ll[4]) == 2:
                
                flag = checkFileType('\t'.join(ll))
                if flag == False:
                    sys.exit()

                llFreq = listFreq(ll[9:])
                hetNum = 0
                
                #缺失基因型仅支持使用"."或者"?"表示，用其他的标记会导致最后的结果出现错误
                missGT = '.' + flag + '.' if '.' in line else '?' + flag + '?'
                for k in llFreq.keys():
                    if k.split(flag)[0] != k.split(flag)[1]:
                        hetNum += llFreq[k]                
                if missGT in llFreq:
                    t = sum(llFreq.values()) - llFreq[missGT]
                    het = hetNum/t if t != 0 else 1                     
                else:
                    het = hetNum/sum(llFreq.values()) if sum(llFreq.values()) != 0 else 0.0
                a0 = '\t'.join(ll[9:]).count("0")
                a1 = '\t'.join(ll[9:]).count("1")
                maf = min(a0,a1)/(a0+a1) if a0+a1>0 else 0.0           
                miss = llFreq[missGT]/sum(llFreq.values()) if missGT in llFreq else 0.0
                if writeDPFile != None or writeReportFile != None:
                    DP = np.sum(DPs)/len(DPs) if len(DPs) != 0 else 1
                if writeHetFile != None:
                    writeHet.write('_'.join(ll[:2]) + '\t' + str(het) + '\n')
                if writeMafFile != None:
                    writeMaf.write('_'.join(ll[:2]) + '\t' + str(maf) + '\n')
                if writeReportFile != None:
                    writeReport.write('\t'.join(['_'.join(ll[:2])+str(maf),\
                     str(het), str(miss), str(DP)]) + '\n')
                if writeDPFile != None:
                    writeDP.write('_'.join(ll[:2]) + '\t' + str(DP) + "\n")
                if het <= maxHetThreshold and het >= minHetThreshold and \
                    maf <= maxMAFThreshold and maf >= minMAFThreshold and miss <= maxMissThreshold:
                        if writeOutFile != None:
                            writeOut.write('\t'.join(ll)+'\n')

    if readInFile != None: readIn.close()
    if writeOutFile != None: writeOut.close()
    if writeHetFile != None: writeHet.close()
    if writeMafFile != None: writeMaf.close()
    if writeReportFile != None: writeReport.close()
    if writeDPFile != None: writeDP.close()

                
if __name__ == "__main__":
    
    from optparse import OptionParser
    optparser = OptionParser(description="A script used to process VCF file!")
    optparser.add_option("-i", "--input", default=None, dest="input", help="input, *.vcf or *.vcf.gz")
    optparser.add_option("-o", "--output", default=None, dest="output", help="output, VCF")
    optparser.add_option("-H", "--het", default=None, dest="het", help="A file you can write the SNPs' het into.")
    optparser.add_option("-m", "--max-miss", default="1", dest="miss", help="max miss threshold, default is 1")
    optparser.add_option("-I", "--individuals", default=None, dest="individual", help="A file contains individuals line by line")
    optparser.add_option("-e", "--min-het", default="0", dest="minHetThreshold", help="min het threshold, default is 0.")
    optparser.add_option("-t", "--max-het", default="1", dest="maxHetThreshold", help="max het threshold, default is 1.")
    optparser.add_option("-a", "--min-maf", default="0", dest="minMAFThreshold", help="min MAF threshold, default is 0.")
    optparser.add_option("-x", "--max-maf", default="1", dest="maxMAFThreshold", help="max MAF threshold, default is 1.")
    optparser.add_option("-f", "--maf-file", default=None, dest="maf", help="maf file.")
    optparser.add_option("-r", "--report-file", default=None, dest="report", help="report file.")
    optparser.add_option("-d", "--DP-file", default=None, dest="DP", help="DP")
    optparser.add_option("-c", "--change-Chromo", default=False, dest="changeChromo", help="True or False. Convert the \"chrA01\" to \"1\".")

    
    (options, args) = optparser.parse_args()
    
    assert options.input != None
    preData(options.input, minHetThreshold=options.minHetThreshold, maxHetThreshold = options.maxHetThreshold,\
            minMAFThreshold=options.minMAFThreshold, maxMAFThreshold=options.maxMAFThreshold, maxMissThreshold=options.miss,\
            writeOutFile=options.output, writeHetFile=options.het, filterIndividualFile=options.individual, writeMafFile=options.maf,\
            writeReportFile=options.report, writeDPFile=options.DP, convert=options.changeChromo)
