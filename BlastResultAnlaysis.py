import pandas as pd
import pyfaidx
import math
import csv
import numpy as np
import sys

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# procedure step 1: obtain csv
fields = ['Chr', 'PamNum', 'OptimalNum', 'gRNA', 'gRNA_Pos']
df = pd.read_csv('Output.csv', skipinitialspace=True, usecols=fields)

outputFileName = "BlastOutput.csv"

rowMaxNum = df.shape[0]     # get the row num of data
gRNA_Str = ""
gRNA_Pos_Str = ""
PamNum = 0
optimalNum = 0
rownum = 0

gRNA_List = []          # stored all founded PAM sequences
gRNA_Pos_List = []      # stored all founded PAM sequence positions
hitNum = 0              # # of hit
hitScore = 0            # the score of hit
seqLen = 0              #  the length of aligned sequence
similarity_percent = 0  # percent of similarity = hitSore / 2* seqLen

BlastMatrix = []

while (rownum <= rowMaxNum - 1):    # extract all sequences to be searched from output file
    Chr = df.loc[rownum, 'Chr']

    if not math.isnan(Chr):  # check not nan string
        PamNum = int(df.loc[rownum, 'PamNum'])
        optimalNum = int(df.loc[rownum, 'OptimalNum'])
        gRNA_Str = str(df.loc[rownum, 'gRNA'])
        gRNA_Pos_Str = str(df.loc[rownum, 'gRNA_Pos'])

        if (gRNA_Str != ""):   # not an empty string
            gRNA_StrList = gRNA_Str.split()

            for gRNA in gRNA_StrList:
                gRNA_Seq = Seq(gRNA)
                if(str(gRNA_Seq) != "nan"):
                    print ("Sequence to be searched: " + str(rownum) + "\t" + str(gRNA_Seq))
                    gRNA_List.append(str(gRNA_Seq))

        if (gRNA_Str != ""):  # not an empty string
            gRNA_Pos_StrList = gRNA_Pos_Str.split()

            for gRNA_Pos in gRNA_Pos_StrList:
                if (str(gRNA_Pos) != "nan"):
                    gRNA_Pos_List.append(str(gRNA_Pos))

    rownum += 1

length = len(gRNA_List)
for i in range(length):
    gRNA_Str = gRNA_List[i]
    gRNA_Pos_Str = gRNA_Pos_List[i]

    fileName = gRNA_Str + "_blast.xml"
    result_handle = open(fileName)
    blast_record = NCBIXML.read(result_handle)

    E_VALUE_THRESH = 4.0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("e value:", hsp.expect)
            if hsp.expect < E_VALUE_THRESH:
                print("\nAlignment for " + gRNA_Str)
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("# of Hits:", alignment.hit_id)
                print("score:", hsp.score)
                print("gaps:", hsp.gaps)
                print("e value:", hsp.expect)

                hitNum = hsp.num_alignments  # # of hit
                hitScore = hsp.score  # the score of hit
                seqLen = len(gRNA_Str)  # the length of aligned sequence
                similarity_percent = int(hitScore*100/(2*seqLen))  # percent of similarity = hitSore / 2* seqLen
                BlastMatrix.append([gRNA_Pos_Str, gRNA_Str, str(hitNum), str(similarity_percent)])

df = pd.DataFrame(BlastMatrix)  # create the dataframe from the result matrix
df.columns = ['gRNA Position', 'gRNA Sequence', '# of Hits', 'Similarity']  # set the name fields

# write to CVS file
df.to_csv(outputFileName)

print("\nResults were written to " + outputFileName)