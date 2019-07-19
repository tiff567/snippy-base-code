import pandas as pd
import pyfaidx
import math
import csv
import numpy as np
import sys

from Bio.Seq import Seq

# import all functions of sequence extraction based on Cas9 type
from PamSeqExtraction import *

# get input from the user
inputFileName = 'Input.csv'
outputFileName = 'output.csv'
Cas9Type = 'NGG'   # supported Cas9 Type: N-GG, N-GRRT/N-GRRN, N-NNNGATT, N-NNNRYAC, N-NNAGAAW
isOptimal = True

shiftPos = 14  # the shift position for locating the PAM
seqLen = 20   # The length of Pam Sequences

# procedure step 1: obtain csv
fields = ['Chr', 'Pos', 'Ref', 'Alt']
df = pd.read_csv(inputFileName, skipinitialspace=True, usecols=fields)

# check the file whether has correct value
rowMaxNum = df.shape[0]     # get the row num of data
rownum = 0
fileEmpty = False

while (rownum <= rowMaxNum - 1):
    Chr = df.loc[rownum, 'Chr']

    if not math.isnan(Chr):  # check not nan string
        fileEmpty = True

        if (int(Chr) < 1) or (int(Chr) > 23):    # Chr is out of range
            sys.exit(" Wrong data file inputted")

    rownum += 1

if (fileEmpty == False):
    sys.exit(" Empty data file inputted")   # empty file

# procedure step 2: import human genome sequences
from pyfaidx import Fasta
hg19 = Fasta('GRCh37.p13.genome.fa')

# procedure step 3: sanity check if there is A or not at the position, if false say something, make table, append, make a false
rownum = 0
sanity = 'YES'
rowMaxNum = df.shape[0]     # get the row num of data

resultMatrix = []  # define a matrix for output results

gRNA_Str = ""       # contain PAM sequences
gRNA_Pos_Str = ""     # contain PAM sequence positioins
PamNum = 0          # contain the # of Pams
optimalNum = 0      # contain the # of optimal Pams

while (rownum <= rowMaxNum - 1):
    Chr = df.loc[rownum, 'Chr']

    if not math.isnan(Chr):  # check not nan string
        Chr = int(df.loc[rownum, 'Chr'])
        Pos = df.loc[rownum, 'Pos']
        Ref = df.loc[rownum, 'Ref']
        Alt = df.loc[rownum, 'Alt']

        # link https://proquest.com/openview/66b90475c5bd328de0363512c8a6ef8c/1?pq-origsite=gscholar&cbl=2045933
        FasRef = hg19['chr'+str(Chr)][Pos:Pos+1]

        # step 1: check whether it is sanity or not
        if str(FasRef) != Ref:    # found unmatched Ref
            sanity = 'NO'
        else:
            sanity = 'YES'

        # step 2: check mutation type A->G or T->C or No
        if Ref == "A" and Alt == "G":   # A->G mutation
            mutation = "AG"
        elif Ref == "T" and Alt == "C":  # T->C mutation
            mutation = "TC"
        else:
            mutation = "noMutation"

        geneSense = "None"

        # step 3 check PAM, optimal and get the gRNA sequence for A->G
        if (mutation == "AG") & (sanity == "YES"):  # AG mutation

            # Based on given Cas9 Type, extract corresponding pam sequences
            # supported Cas9 Type: N-GG, N-GRRT/N-GRRN, N-NNNGATT, N-NNNRYAC, N-NNAGAAW
            if Cas9Type == 'NGG' or Cas9Type.lower() == 'spcas9':
                myResult = get_pam_sequence_ag_spcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif (Cas9Type =='NGRRT') or (Cas9Type == 'NGRRN') or (Cas9Type.lower() == 'sacas9'):
                myResult = get_pam_sequence_ag_sacas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type =='NNNNGATT' or Cas9Type.lower() == 'nmecas9':
                myResult = get_pam_sequence_ag_nmecas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type =='NNNNRYAC' or Cas9Type.lower() == 'cjcas9':
                myResult = get_pam_sequence_ag_cjcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type =='NNNAGAAW' or Cas9Type.lower() == 'stcas9':
                myResult = get_pam_sequence_ag_stcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            else:
                sys.exit("Unsupported CAS9 Type. Please try other type again !!")

            gRNA_Str, gRNA_Pos_Str, PamNum, optimalNum = myResult

            if isOptimal:
                if(optimalNum > 0):
                    geneSense = "Sense"
                    print("Sense gRNA sequences:")
                    print("row #:", rownum, gRNA_Str)

            # write done A->G mutation information
            resultMatrix.append([sanity, mutation, str(PamNum), str(optimalNum), gRNA_Str, gRNA_Pos_Str, geneSense])

         # step 4 check PAM, optimal and get the gRNA sequence for T->C

        elif (mutation == "TC") & (sanity == "YES"):  # T->C mutation

            # Based on given Cas9 Type, extract corresponding pam sequences
            # supported Cas9 Type: N-GG, N-GRRT/N-GRRN, N-NNNGATT, N-NNNRYAC, N-NNAGAAW
            if Cas9Type == 'NGG' or Cas9Type.lower() == 'spcas9':
                myResult = get_pam_sequence_tc_spcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif (Cas9Type == 'NGRRT') or (Cas9Type == 'NGRRN') or (Cas9Type.lower() == 'sacas9'):
                myResult = get_pam_sequence_tc_sacas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type == 'NNNNGATT' or Cas9Type.lower() == 'nmecas9':
                myResult = get_pam_sequence_tc_nmecas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type == 'NNNNRYAC' or Cas9Type.lower() == 'cjcas9':
                myResult = get_pam_sequence_tc_cjcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            elif Cas9Type == 'NNNAGAAW' or Cas9Type.lower() == 'stcas9':
                myResult = get_pam_sequence_tc_stcas9(hg19, Chr, Pos, shiftPos, seqLen, isOptimal)
            else:
                sys.exit("Unsupported CAS9 Type. Please try other type again !!")

            gRNA_Str, gRNA_Pos_Str, PamNum, optimalNum = myResult

            if isOptimal:
                if (optimalNum > 0):
                    geneSense = "Sense"
                    print("Sense gRNA sequences:")
                    print("row #:", rownum, gRNA_Str)

            # write done T->C mutation information
            resultMatrix.append([sanity, mutation, str(PamNum), str(optimalNum), gRNA_Str, gRNA_Pos_Str, geneSense])

        else:
             # write done none mutation information
            geneSense = "NaN"
            PamNum = 0
            optimalNum = 0
            gRNAStr = ''

            resultMatrix.append([sanity, mutation, str(PamNum), str(optimalNum), gRNAStr, geneSense])
    rownum += 1


df1 = pd.read_csv('Input.csv', skipinitialspace=True)   # get the original dataframe

df2 = pd.DataFrame(resultMatrix)  # create the dataframe from the result matrix
df2.columns = ['IsSanity', 'MutationType', 'PamNum', 'OptimalNum', 'gRNA', 'gRNA_Pos','SenseType']  # set the name fields

dfnew = pd.concat([df1, df2], axis=1)    # concatenate the two dataframes by columns

# write to CVS file
dfnew.to_csv(outputFileName)

print("\nResults were written to output.csv!")


