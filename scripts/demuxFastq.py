import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gzip
#2442687652
libList = pd.read_csv('data/AggrFileListcDNAwithLib.csv')
goodCells = pd.read_csv('data/ClassifiedDemuxedBarcodes_July23_2020.csv.gz')
goodCells=goodCells.loc[(goodCells.barcodes != 'Negative') & (goodCells.barcodes != 'Doublet')]
goodCells.reset_index(drop=True)

## for each fastq library, iterate and make good cells
goodCellSeries=goodCells.cellIds.str.split('.').to_list()
goodCellSeries
cellMaster=defaultdict(list)
for i, j, in goodCellSeries:
    cellMaster[j].append(i)


def gzipHandle(fileName):
    if '.gz' in fileName:
        fileOut = gzip.open(fileName, 'rt')
    else:
        fileOut = open(fileName, 'rt')
    return fileOut

def fileOutOp(sampleBC):
    R1out = os.path.join('outs', sampleBC+'_'+'R1.fastq')
    R2out = os.path.join('outs', sampleBC+'_'+'R2.fastq')
    if os.path.exists(R1out) and os.path.exists(R2out):
        R1outHandle = open(R1out, 'a')
        R2outHandle = open(R2out, 'a')
    else:
        R1outHandle = open(R1out, 'w')
        R2outHandle = open(R2out, 'w')
    return [R1outHandle, R2outHandle]


def processFastq(R1, R2, libN, goodCells, cellMaster):
    #libN = '1'
    R1File = gzipHandle(R1)
    R2File = gzipHandle(R2)
    r1String = ''
    r2String = ''
    umiMaster = defaultdict(str)
    writeState = False
    fastqRec = 0
    for lineN, (r1, r2) in enumerate(zip(R1File, R2File)):
        if lineN % 10000 == 0:
            print('{} PROCESSED LINES {} AND {} FASTQ Recs'.format(R1, lineN, fastqRec))
        if lineN % 4 == 1: # every 2nd lines
            fastaLineR1 = r1.strip() ## Cell barcode will appear at 0:16
            #fastaLineR2 = r2.strip()
            cell, umiInd= [0, 16], [16, 28]
            cell_start, cell_end = cell
            umi_start, umi_end = umiInd
            cellBC, umi = fastaLineR1[cell_start:cell_end], fastaLineR1[umi_start:umi_end]
            if umi in umiMaster:
                writeState = True
                sampleBC = umiMaster.get(umi)
                if sampleBC:
                    R1out, R2out = fileOutOp(sampleBC)
                    r1String += r1
                    r2String += r2  
            else:
                if cellBC in cellMaster.get(libN):
                    writeState = True
                    sampleBC=goodCells.loc[goodCells.cellIds == cellBC+'.'+libN, 'barcodes'].to_list()[0]
                    if sampleBC:# == "ACATGCGT":
                        umiMaster[umi] = sampleBC
                        #assert len(sampleBC) == 1
                        ## make file handler function here
                        R1out, R2out = fileOutOp(sampleBC)
                        r1String += r1
                        r2String += r2  
                else:
                    #r1String = r2String = ''
                    writeState = False
        elif lineN % 4 == 3: # end fast rec
            if writeState is True:
                #print(writeState)
                r1String += r1
                r2String += r2  
                R1out.write(r1String)
                R2out.write(r2String)
                r1String = r2String = '' # reset the string
                writeState = False
                fastqRec += 1
            else:
                r1String = r2String = ''
        else:
            r1String += r1
            r2String += r2     
    R1File.close()
    R2File.close()
    R1out.close()
    R2out.close()

R1 = 'data/cDNA8_S88_L004_R1_001.fastq.gz'
R2 = 'data/cDNA8_S88_L004_R2_001.fastq.gz'
libN = '16'

processFastq(R1, R2, libN, goodCells, cellMaster)